/*
 *  Main author:
 *     Mikael Östlund <oestlund.mikael@gmail.com>
 *
 *  Copyright:
 *     Mikael Östlund, 2017
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <gecode/int.hh>

using namespace Gecode;
using namespace Gecode::Int;

/* The FR class represents forbidden regions using two arrays, min and max */
class FR {
public:
  int *min; // origin per dimension
  int *max; // maximal per dimension
};

/* The Object class represents a hyperrectangle object */
class Object {
public:
  ViewArray<IntView> x; // Object origins
  int *l; // Object lengths
  int id; // Object id
  int *support_min; // 1D matrix representation for keeping track of supported points for pruneMin
  int *support_max; // 1D matrix representation for keeping track of supported points for pruneMax
  
  // Arrays used in ADVISOR optimisation
  int *rfrb;
  int *rfre;
  int *d_size;

  // Boolean used in SEPARATE optimisation
  bool fixed;
  
  bool isSame(Object *);
};

bool Object::isSame(Object *o) {
  return id == o->id;
}

/* The ForbiddenRegions class is used for storing forbidden regions and retrieving them using RR scheduling */
class ForbiddenRegions {
private:
  int RRpos; // Current RR position in GetFR
  int length; // Amount of FRs in collection
public:
  Support::DynamicArray<FR*, Region> collection; // TODO: this will be private in the future
public:
  int size();
  void insert(FR*);
  FR* getRR();
  void resetRR();
  ForbiddenRegions(Region*);
};


ForbiddenRegions::ForbiddenRegions(Region *r) : collection(*r), length(0), RRpos(0) {}

FR* ForbiddenRegions::getRR() {
  return collection[(RRpos++ % length)];
}

void ForbiddenRegions::resetRR() {
  RRpos = 0;
}

int ForbiddenRegions::size(void) {
  return length;
}

void ForbiddenRegions::insert(FR *f) {
  collection[length] = f;
  ++length;
}

/* The OBJECTS class, used for storing Objects */
class OBJECTS {
private:
  int length;
public:
  Support::DynamicArray<Object*, Space> collection;
  int size();
  void insert(Object*);
  void prettyPrint();

  OBJECTS(Space& s) : collection(s), length(0) {};
};

void OBJECTS::prettyPrint() {
  for (int i = 0; i < size(); i++) {
    Object *o = collection[i];
    std::cout << o->x << "\n";
  }
}

int OBJECTS::size(void) {
  return length;
}

void OBJECTS::insert(Object *f) {
  collection[length] = f;
  ++length;
}

// The diffn propagator
class Diffn : public Propagator {
protected:
  OBJECTS *Objects; // Objects being filtered
  int dimensions; // Number of dimensions of the problem

  int *maxl;

  class ViewAdvisor : public Advisor {
  public:
    //Int::IntView x;
    int dim; // What dimension in the corresponding viewarray is this advisor responsible for?
    int i; // The corresponding object
    ViewAdvisor(Space& home, Propagator& p, 
                Council<ViewAdvisor>& c, int d, int idx) 
      : Advisor(home,p,c), dim(d), i(idx) {
    }
    ViewAdvisor(Space& home, bool share, ViewAdvisor& a)
      : Advisor(home,share,a), dim(a.dim), i(a.i) {
    }
    void dispose(Space& home, Council<ViewAdvisor>& c) {
      Advisor::dispose(home,c);
    }
  };
  
  Council<ViewAdvisor> c;
  
  // True iff f and o overlap in dimension d
  bool overlaps(int min, int max, Object *o, int d) {
    if ((o->x[d].max() < min) || (o->x[d].min() > max)) {
      return false;
    }
    return true;
  }

  // Generates relative FRs to the object o and stores them in F
  void genOutBoxes(Region &r, ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];

      if (!other->isSame(o)) { // For every object <> o
        FR *f = (FR *) r.alloc<FR>(1);
        f->min = r.alloc<int>(k*2); // Making sure f->max also fits
        f->max = &(f->min[k]);
        bool exists = true; // Assume f exists

        for (int d = 0; d < k; d++) { // For every dimension d
          const int min = other->x[d].max() - o->l[d] + 1;
          const int max = other->x[d].min() + other->l[d] - 1;
          if ((min <= max) && overlaps(min, max, o, d)) {
            f->min[d] = min;
            f->max[d] = max;
          } else {
            exists = false;
            break; // Break loop, other is not interfering with o
          }
        }

        // Only add forbidden region if it exists
        if (exists) {
          F->insert(f);
        }
      }
    }
  }

  // Checks if object o can be skipped during relative FR generation
  bool canSkip(Object *o, int k) {
    for (int i = 0; i < k; i++) {
      if (o->d_size[i] > (maxl[i] - 2)) {
        return true;
      }
    }
    return false;
  }

  /* Tries to coalesce f0 and f1.
   * Return cases:
   * 0: coalescing not possible
   * 1: f0 includes f1 (completely overlaps)
   * 2: f1 includes f2 (completely overlaps)
   * 3: one of the two can be extended in one dimension to represent both FRs,
   *    f0 is changed in place to represent both.
   */
  int coalesce(FR *f0, FR *f1, int k) {
    // TODO: introduce enumerator
    int trend = 0; // What case did the previous dimensions fullfil?
    int e = 0; // In case 3, what dimension can be extended?

    for (int d = 0; d < k; d++) {
      if (f0->max[d] + 1 < f1->min[d] || f0->min[d] > f1->max[d] + 1) { /* Do not coalesce */
        return 0;
      } else if (f0->min[d] == f1->min[d] && f0->max[d] == f1->max[d]) { /* Equal in dimension d */
        // no-op
      } else if (f0->min[d] <= f1->min[d] && f0->max[d] >= f1->max[d]) { /* f1 \subset f0*/
        if (trend == 0 || trend == 1) {
          trend = 1;
        } else {
          return 0; // The FRs cannot be coalesced
        }
      } else if (f0->min[d] >= f1->min[d] && f0->max[d] <= f1->max[d]) { /* f0 \subset f1*/
        if (trend == 0 || trend == 2) {
          trend = 2;
        } else {
          return 0; // The FRs cannot be coalesced
        }
      } else { /* The FRs intersect eachother in dimension d */
        e = d;
        if (trend == 0) {
          trend = 3;
        } else {
          return 0; // The FRs cannot be coalesced
        }
      }
    }

    switch (trend) {
    case 0: /* They are equal*/
      return 1; // Let f0 represent them both
    case 1: /* f0 includes f1 */
      return 1;
    case 2: /* f1 includes f0 */
      return 2;
    case 3: /* They overlap in dimension e, let f0 represent them both */
      f0->min[e] = std::min(f0->min[e], f1->min[e]);
      f0->max[e] = std::max(f0->max[e], f1->max[e]);
      return 1;
    }
  }

  /* Generates forbidden regions for object o given objects in O, merges forbidden regions where possible */
  void genOutBoxesMerge(Region& r, ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    Support::DynamicQueue<FR*, Region> Q(r); // Queue for temporary storage of merge candidates
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];
      
      if (!other->isSame(o)) { // For every other object <> o
        if (canSkip(other, k)) {
          //std::cout << "Skipping " << other->id << "\n";
          continue;
        }

        FR *f = (FR *) r.alloc<FR>(1);
        f->min = r.alloc<int>(k*2);
        f->max = &(f->min[k]);

        bool exists = true; // Assume f exists

        for (int d = 0; d < k; d++) { // For every dimension d
          const int min = other->rfrb[d] - o->l[d];//other->x[d].max() - o->l[d] + 1;
          const int max = other->rfre[d];//other->x[d].min() + other->l[d] - 1;
          if ((min <= max) && overlaps(min, max, o, d)) {
            f->min[d] = min;
            f->max[d] = max;
          } else {
            exists = false;
            break; // Break loop, other is not interfering with o
          }
        }

        if (exists) {
          while (exists && !Q.empty()) {
            FR *f0 = Q.pop();
            switch (coalesce(f0, f, k)) { /* Try to coalesce the FRs */
            case 0:  /* No coalescing possible */
              Q.push(f0);
              exists = false;
              break;
            case 1: /* f0 subsumes f */
              f = f0;
              break;
            case 2: /* f subsumes f0 */
              break;
            }
          }
          Q.push(f);
        }
      }
    }

    while (!Q.empty()) {
      F->insert(Q.pop());
    }
  }

  // Returns true iff c is feasible according to outbox
  bool isFeasible(FR *outbox, int k, int *c) {
    for (int d = 0; d < k; d++) {
      if ((c[d] < outbox->min[d]) || (c[d] > outbox->max[d])) {
        return true;
      }
    }
    return false;
  }

  // Returns next infeasible FR if one exists, NULL otherwise
  FR *getFR(int k, int *c, ForbiddenRegions *actrs) {
    for (int i = 0; i < actrs->size(); i++) {
      FR *f = actrs->getRR(); // Using RR for getting the next FR
      if (!isFeasible(f, k, c)) { // forbidden region for c found
        return f; // f is an infeasible FR
      }
    }
    return NULL; // return NULL if no infeasible FR was found
  }

  // Function for pruning the lower bound of object o in dimension d
  ExecStatus pruneMin(Home home, Object *o, int d, int k, ForbiddenRegions *F) {
    // SUPPORT optimisation
    bool supported = true;

    /* SUPPORT optimisation */
    for (int j = 0; j < k; j++) {
      if (!o->x[j].in(o->support_min[d*k+j])) {
        supported = false;
        break;
      }
    }

    /* SUPPORT optimisation */
    // Considering the d:th row in the support_min "matrix"
    if (supported && !getFR(k, &(o->support_min[d*k]), F)) {
      return ES_FIX;
    }

    bool b = true;    
    Region r(home); // Region for storage

    /* Init c and n*/
    int *c = r.alloc<int>(k); // sweep-point
    // init c
    for (int i = 0; i < k; i++) {
      c[i] = o->x[i].min();
    }

    int *n = r.alloc<int>(k); // jump vector
    // init vector
    for (int i = 0; i < k; i++) {
      n[i] = o->x[i].max() + 1;
    }

    // Get next FR
    FR *currentF = getFR(k, c, F);

    // While we have not failed and c is infeasible
    while (b and currentF) {
      // update jump-vector
      for (int i = 0; i < k; i++) { // TODO: abstract to updatevector procedure
        n[i] = std::min(n[i], currentF->max[i] + 1);
      }
      
      b = false; // No new point to jump to yet (assume failure)
      
      // TODO: Abstract to adjust procedure
      // Adjust the sweep-point based on jump vector information
      for (int j = k - 1; j >= 0; j--) {
        int r = (j + d) % k;         // Consider rotation
        c[r] = n[r];                 // Use n vector to jump
        n[r] = o->x[r].max() + 1;    // Reset component of n after jumping
        if (c[r] <= o->x[r].max()) { // Jump target found?
          b = true;                  // target found
          break;
        } else {
          c[r] = o->x[r].min();      // reset component of c, dimension r exhausted
        }
      }
      // Update currentF and check if it is infeasible still
      currentF = getFR(k, c, F);
    }

    if (b) { 
      /* SUPPORT optimisation */
      for (int j = 0; j < k; j++) {
        o->support_min[d*k+j] = c[j];
      }

      ModEvent me = o->x[d].gq(home, c[d]); // prune o
      if (me_failed(me)) {
        return ES_FAILED;
      } else if (me_modified(me)) {
        o->support_min[d * k + d] = o->x[d].min(); // In case c[d] was in a hole
        o->rfre[d] = o->x[d].min() + o->l[d] - 1;
        o->d_size[d] = o->x[d].max() - o->x[d].min() - o->l[d];
        return ES_NOFIX;
      } else {
        return ES_FIX;
      }
    } else {
      return ES_FAILED;
    }
  }

  // Function for pruning the upper bound of object o in dimension d
  ExecStatus pruneMax(Home home, Object *o, int d, int k, ForbiddenRegions *F) {
    // SUPPORT optimisation
    bool supported = true;

    /* SUPPORT optimisation */
    for (int j = 0; j < k; j++) {
      if (!o->x[j].in(o->support_max[d*k+j])) {
        supported = false;
        break;
      }
    }

    /* SUPPORT optimisation */
    // Considering the d:th row in the support_max "matrix"
    if (supported && !getFR(k, &(o->support_max[d*k]), F)) {
      return ES_FIX;
    }

    bool b = true;
    Region r(home);

    /* Init c and n*/
    int *c = r.alloc<int>(k); // sweep-point
    // init c
    for (int i = 0; i < k; i++) {
      c[i] = o->x[i].max();
    }    

    int *n = r.alloc<int>(k); // jump-vector
    // initvector
    for (int i = 0; i < k; i++) {
      n[i] = o->x[i].min() - 1;
    }

    // Get next FR
    FR *currentF = getFR(k, c, F);

    // While we have not failed and c is infeasible
    while (b and currentF) {
      for (int i = 0; i < k; i++) { // TODO: abstract to updatevector procedure
        n[i] = std::max(n[i], currentF->min[i] - 1);
      }

      b = false; // No new point to jump to yet (assume failure)
      
      // TODO: Abstract to adjust procedure
      for (int j = k - 1; j >= 0; j--) {
        int r = (j + d) % k;
        c[r] = n[r];              // Use n vector to jump
        n[r] = o->x[r].min() - 1;    // Reset component of n after jumping
        if (c[r] >= o->x[r].min()) { // Jump target found?
          b = true;               // target found
          break;
        } else {
          c[r] = o->x[r].max();      // reset component of c, dimension r exhausted
        }
      }

      // Update currentF and check if it is infeasible still
      currentF = getFR(k, c, F);
    }

    if (b) { 
      /* SUPPORT optimisation */
      for (int j = 0; j < k; j++) {
        o->support_max[d*k+j] = c[j];
      }

      ModEvent me = o->x[d].lq(home, c[d]); // prune o
      if (me_failed(me)) {
        return ES_FAILED;
      } else if (me_modified(me)) {
        o->support_max[d*k+d] = o->x[d].max(); // In case c[d] was in a hole
        o->rfrb[d] = o->x[d].max() + 1;
        o->d_size[d] = o->x[d].max() - o->x[d].min() - o->l[d];
        return ES_NOFIX;
      } else {
        return ES_FIX;
      }
    } else {
      return ES_FAILED; // No feasible sweep point found, report failure
    } 
  }

  // R is a collection of rectangles participating in the problem
  ExecStatus filter(Space &home, int k) {
    bool nonfix = true;
    bool allfixed = true; // Used for detecting subsumption

    while (nonfix) {
      nonfix = false;
      allfixed = true;
      for (int i = 0; i < Objects->size(); i++) {
        Object *o = Objects->collection[i];

        // SEPARATE: if o is fixed and checked, then skip it.
        if (o->fixed) {
          continue; 
        }

        Region r(home); // TODO: one region per object or one for all objects? ("keep scope small")
        ForbiddenRegions F(&r); 
        genOutBoxesMerge(r, &F, Objects, k, o);

        if (o->x.assigned()) {
          int c[k];
          for (int i = 0; i < k; i++) {
            c[i] = o->x[i].val();
          }
          if (getFR(k, c, &F)) {
            return ES_FAILED;
          } else {
            F.resetRR();
          }
        } else {
          for (int d = 0; d < k; d++) {
            ExecStatus pMinStatus = ES_FIX;
            if (!o->x[d].assigned()) {
              pMinStatus = pruneMin(home, o, d, k, &F);
              if (pMinStatus == ES_FAILED) { // No feasible minimal point found, report failure
                return ES_FAILED;
              }
            }
            F.resetRR(); // Reset RR position
            ExecStatus pMaxStatus = ES_FIX;
            if (!o->x[d].assigned()) {
              pMaxStatus = pruneMax(home, o, d, k, &F);
              if (pMaxStatus == ES_FAILED) { // No feasible maximal point found, report failure
                return ES_FAILED;
              }
            }
            if (pMinStatus == ES_NOFIX || pMaxStatus == ES_NOFIX) { // We pruned a bound, not at fixpoint
              nonfix = true;
            }
          }
        }

        if (!o->x.assigned()) {
          allfixed = false;
        } else {
          o->fixed = true;
        }
      }
    }

    // If all objects are fixed and we have not failed, we can subsume
    if (allfixed) {
      return home.ES_SUBSUMED(*this);
    }

    return ES_FIX; // Filter is a fixpoint loop, thus report at fixpoint
  }

public:
  // Create propagator and initialize
  Diffn(Home home, // Constructor for 2D
        ViewArray<IntView>& x0,int w0[],
        ViewArray<IntView>& y0,int h0[])
    : Propagator(home), c(home)
  {
    Objects = (OBJECTS*)((Space &) home).ralloc(sizeof(OBJECTS));
    new(Objects) OBJECTS((Space &) home);

    maxl = ((Space&) home).alloc<int>(2);
    maxl[0] = -1;
    maxl[1] = -1;

    // Create corresponding objects for the ViewArrays and arrays
    for (int i = 0; i < x0.size(); i++) {
      Object *o = ((Space&) home).alloc<Object>(1);
      o->fixed = false;

      o->l = ((Space&) home).alloc<int>(2);
      o->support_min = (int *)((Space&) home).ralloc(sizeof(int) * 2 * 2);
      o->support_max = (int *)((Space&) home).ralloc(sizeof(int) * 2 * 2);
      o->x = ViewArray<IntView>((Space&) home, 2);

      o->x[0] = x0[i];
      o->x[1] = y0[i];
      o->l[0] = w0[i];
      o->l[1] = h0[i];
      o->id = i;

      maxl[0] = std::max(maxl[0], o->l[0]);
      maxl[1] = std::max(maxl[0], o->l[1]);

      o->support_min[0] = o->x[0].min();
      o->support_min[1] = o->x[1].min();
      o->support_min[2] = o->x[0].min();
      o->support_min[3] = o->x[1].min();

      o->support_max[0] = o->x[0].max();
      o->support_max[1] = o->x[1].max();
      o->support_max[2] = o->x[0].max();
      o->support_max[3] = o->x[1].max();

      o->rfre = ((Space&) home).alloc<int>(2*2);
      o->rfrb = &(o->rfre[2]);

      o->rfrb[0] = o->x[0].max() + 1;
      o->rfre[0] = o->x[0].min() + o->l[0] - 1;
      o->rfrb[1] = o->x[1].max() + 1;
      o->rfre[1] = o->x[1].min() + o->l[1] - 1;

      o->d_size = ((Space&) home).alloc<int>(2);

      o->d_size[0] = o->x[0].max() - o->x[0].min() - o->l[0];
      o->d_size[1] = o->x[1].max() - o->x[1].min() - o->l[1];

      o->x[0].subscribe(home,*new (home) ViewAdvisor(home,*this,c,0,i));
      o->x[1].subscribe(home,*new (home) ViewAdvisor(home,*this,c,1,i));

      Objects->insert(o);
    }

    dimensions = 2;

    Int::IntView::schedule(home, *this, Int::ME_INT_BND); // Schedule the propagator
    home.notice(*this, AP_DISPOSE); // Make sure dispose function is called on Space destruction
  }

  // Post no-overlap propagator
  static ExecStatus post(Home home,
                         ViewArray<IntView>& x, int w[],
                         ViewArray<IntView>& y, int h[]) {
    // Only if there is something to propagate
    if (x.size() > 1)
      (void) new (home) Diffn(home,x,w,y,h);
    return ES_OK;
  }

  // Dispose propagator and return its size
  virtual size_t dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    c.dispose(home);

    (void) Propagator::dispose(home);
    return sizeof(*this);
  }
    
  // Copy constructor during cloning
  Diffn(Space& home, bool share, Diffn& p)
    : Propagator(home,share,p) {
    dimensions = p.dimensions;
    Objects = (OBJECTS*)((Space &) home).ralloc(sizeof(OBJECTS));
    new(Objects) OBJECTS((Space &) home);

    maxl = ((Space&) home).alloc<int>(dimensions);
    
    for (int j = 0; j < dimensions; j++) {
      maxl[j] = p.maxl[j];
    }

    for (int i = 0; i < p.Objects->size(); i++) {
      Object *pObj = p.Objects->collection[i];
      Object *o = ((Space&) home).alloc<Object>(1);
      o->fixed = pObj->fixed;

      o->support_min = (int *) home.ralloc(sizeof(int) * dimensions * dimensions); // 1D representation of matrix
      o->support_max = (int *) home.ralloc(sizeof(int) * dimensions * dimensions); // 1D representation of matrix

      o->l = home.alloc<int>(dimensions*4); // Make sure memory block fits 3 more arrays of identical size
      o->rfre = &(o->l[dimensions]);
      o->rfrb = &(o->rfre[dimensions]);
      o->d_size = &(o->rfrb[dimensions]);

      o->x.update(home, share, pObj->x);

      for (int j = 0; j < dimensions; j++) {
        o->l[j] = pObj->l[j];

        for (int d = 0; d < dimensions; d++) {
          o->support_min[j * dimensions + d] = pObj->support_min[j * dimensions + d];
          o->support_max[j * dimensions + d] = pObj->support_max[j * dimensions + d];
        }
        
        o->rfre[j] = pObj->rfre[j];
        o->rfrb[j] = pObj->rfrb[j];

        o->d_size[j] = pObj->d_size[j];
      }
      o->id = pObj->id;
      Objects->insert(o);
    }

    c.update(home, share, p.c);
  }
  // Create copy during cloning
  virtual Propagator* copy(Space& home, bool share) {
    return new (home) Diffn(home,share,*this);
  }
    
  // Return cost (defined as expensive quadratic)
  virtual PropCost cost(const Space&, const ModEventDelta&) const {
    return PropCost::quadratic(PropCost::HI,dimensions*Objects->size());
  }

  virtual void reschedule(Space& home) {
    Int::IntView::schedule(home, *this, Int::ME_INT_BND);
  }

  virtual ExecStatus advise(Space& home, Advisor& a, const Delta& d) {
    int dim = (static_cast<ViewAdvisor&>(a)).dim;
    int i = (static_cast<ViewAdvisor&>(a)).i;
    Object *o = Objects->collection[i];

    // update values since view changed
    o->rfrb[dim] = o->x[dim].max() + 1; 
    o->rfre[dim] = o->x[dim].min() + o->l[dim] - 1;
    o->d_size[dim] = o->x[dim].max() - o->x[dim].min() - o->l[dim];

    return ES_NOFIX;
  }
    
  // Perform propagation
  virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
    return filter((Home) home, dimensions);
  }
 
};

/*
 * Post the constraint that the rectangles defined by the coordinates
 * x and y and width w and height h do not overlap.
 *
 * This is the function that you will call from your model. The best
 * is to paste the entire file into your model.
 */
void nonOverlapping(Home home,
                    const IntVarArgs& x, const IntArgs& w,
                    const IntVarArgs& y, const IntArgs& h) {
  // Check whether the arguments make sense
  if ((x.size() != y.size()) || (x.size() != w.size()) ||
      (y.size() != h.size()))
    throw ArgumentSizeMismatch("nooverlap");
  // Never post a propagator in a failed space
  if (home.failed()) return;
  // Set up array of views for the coordinates
  ViewArray<IntView> vx(home,x);
  ViewArray<IntView> vy(home,y);
  // Set up arrays (allocated in home) for width and height and initialize
  int* wc = static_cast<Space&>(home).alloc<int>(x.size());
  int* hc = static_cast<Space&>(home).alloc<int>(y.size());
  for (int i=x.size(); i--; ) {
    wc[i]=w[i]; hc[i]=h[i];
  }
  // If posting failed, fail space
  if (Diffn::post(home,vx,wc,vy,hc) != ES_OK)
    home.fail();
}


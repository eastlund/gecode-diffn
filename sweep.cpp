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

class FR {
public:
  int *min; // origin per dimension
  int *max; // maximal per dimension
};

/* An object consists of the object origin in each dimension d, the
 * length in dimension d, its id, and whether or not it is fixed,
 * i.e., if it needs to be considered while filtering.
 */
class Object {
public:
  ViewArray<IntView> x; // Object origins
  int *l; // Object lengths
  int id;
  int **support_min; // Arrays keeping track of supported points for pruneMin
  int **support_max; // Arrays keeping track of supported points for pruneMax

  int *rfrb;
  int *rfre;

  int *d_size;
  
  bool isSame(Object *);

  //Object() : fixed(false) {}
};

bool Object::isSame(Object *o) {
  return id == o->id;
}

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
class NonOverlapping : public Propagator {
protected:
  OBJECTS *Objects; // Objects being filtered
  int dimensions; // Number of dimensions of the problem

  int pivot;

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
  
  // True if f and o overlap in dimension d
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
        f->min = r.alloc<int>(k);
        f->max = r.alloc<int>(k);
        //f->max = home.realloc<int>(f->min, dimensions, dimensions*2); // This could give better cache util but seems to give worse performance
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

  // 0 if no merge possible
  // 1 if B can be removed
  // 2 if A can be removed
  int completelyOverlaps(FR *A, FR *B, int k) {
    int r = 0;
    for (int d = 0; d < k; d++) {
      if (A->min[d] <= B->min[d] && A->max[d] >= B->max[d]) {
        if (r == 2) {
          return 0;
        }
        r = 1; // A completely overlaps B in dimension d
      } else if (A->min[d] >= B->min[d] && A->max[d] <= B->max[d]) {
        if (r == 1) {
          return 0;
        }
        r = 2; // B completely overlaps A in dimension d
      } else {
        return 0; // No complete overlap
      }
    }
    return r;
  }

  // Is A and B combined a forbidden region? If so, let A represent them both.
  // Should only be called if completelyOverlaps can't merge regions (as combine makes assumptions based on this)
  bool combine(FR *A, FR *B, int k) { // TODO: check if this actually ever combine something
    int differentSizes = 0;
    int differentDim = 0;
    for (int d = 0; d < k; d++) {
      if (not ((A->min[d] == B->min[d]) && (A->max[d] == B->max[d]))) {
        differentSizes++; // Different in current dimension
        differentDim = d; // Different in dimension d
      }
    }

    if (differentSizes > 1) {
      return false; // A and B differ in size in more than one dimension, cannot be merged
    }

    if ((A->min[differentDim] <= B->min[differentDim]) && (A->max[differentDim] >= B->min[differentDim])) {
      // Change A inplace
      A->max[differentDim] = B->max[differentDim]; // B->max must be larger than A->max before this!
      return true;
    } else if ((B->min[differentDim] <= A->min[differentDim]) && (B->max[differentDim] >= A->min[differentDim])) {
      A->min[differentDim] = B->min[differentDim]; // B->min must be smaller than A->min before this!
      return true;
    }
    return false;
  }

  bool canSkip(Object *o, int k) {
    for (int i = 0; i < k; i++) {
      if (o->d_size[i] > (maxl[i] - 2)) {
        return true;
      }
    }
    return false;
  }

  void genOutBoxesMerge(Region& r, ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    int inQueue = 0;
    FR *fst = NULL;
    FR *snd = NULL;
    //std::cout << "--------------------------\n";
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];
      
      if (canSkip(other, k)) {
        //std::cout << "Skipping " << other->id << "\n";
        continue;
      }

      if (!other->isSame(o)) { // For every other object <> o
        FR *f = (FR *) r.alloc<FR>(1);
        f->min = r.alloc<int>(k);
        f->max = r.alloc<int>(k);
        //f->max = home.realloc<int>(f->min, dimensions, dimensions*2); // This could give better cache util

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
          if (inQueue == 2) { // Check if previous 2 FRs can be merged
            int r = completelyOverlaps(fst, snd, k); // TODO: this should be a switch-case

            if (r == 0) { // fst and snd don't completely overlap one another
              if (combine(fst, snd, k)) { // Can they be combined in another way?
                snd = NULL; // snd was merged into fst
                inQueue--;
              } else {
                F->insert(fst); // No merge possible, remove fst from queue
                fst = snd; // Advance queue one step
                snd = NULL;
                inQueue--;
              }
            } else if (r == 1) { // snd is inside fst
              snd = NULL; // Remove snd from queue
              inQueue--;
            } else { // fst is inside snd
              fst = snd; // Advance queue one step
              snd = NULL;
              inQueue--;
            }
          }
          // Add new FR f to queue
          if (fst == NULL) {
            fst = f;
          } else {
            snd = f;
          }
          inQueue++;
        }
      }
    }

    if (fst != NULL) {
      F->insert(fst);
    }

    if (snd != NULL) {
      F->insert(snd);
    }
  }

  // Returns true if c is feasible according to outbox, otherwise false
  bool isFeasible(FR *outbox, int k, int *c) {
    for (int d = 0; d < k; d++) {
      if ((c[d] < outbox->min[d]) || (c[d] > outbox->max[d])) {
        return true;
      }
    }
    return false;
  }

  // Get next infeasible FR
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
    bool supported = true;

    for (int j = 0; j < k; j++) {
      if (!o->x[j].in(o->support_min[d][j])) {
        supported = false;
        break;
      }
    }

    if (supported && !getFR(k, o->support_min[d], F)) {
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
        int r = (j + d) % k;      // Consider rotation
        c[r] = n[r];              // Use n vector to jump
        n[r] = o->x[r].max() + 1;    // Reset component of n after jumping
        if (c[r] <= o->x[r].max()) { // Jump target found?
          b = true;               // target found
          break;
        } else {
          c[r] = o->x[r].min();      // reset component of c, dimension r exhausted
        }
      }
      // Update currentF and check if it is infeasible still
      currentF = getFR(k, c, F);
    }

    if (b) { 
      for (int j = 0; j < k; j++) {
        o->support_min[d][j] = c[j];
      }

      ModEvent me = o->x[d].gq(home, c[d]); // prune o
      if (me_failed(me)) {
        return ES_FAILED;
      } else if (me_modified(me)) {
        o->support_min[d][d] = o->x[d].min(); // In case c[d] was in a hole
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
    bool supported = true;
    for (int j = 0; j < k; j++) {
      if (!o->x[j].in(o->support_max[d][j])) {
        supported = false;
        break;
      }
    }

    if (supported && !getFR(k, o->support_max[d], F)) {
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
      for (int j = 0; j < k; j++) {
        o->support_max[d][j] = c[j];
      }

      ModEvent me = o->x[d].lq(home, c[d]); // prune o
      if (me_failed(me)) {
        return ES_FAILED;
      } else if (me_modified(me)) {
        o->support_max[d][d] = o->x[d].max(); // In case c[d] was in a hole
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
    bool allfixed = true;

    while (nonfix) {
      nonfix = false;
      allfixed = true;
      for (int i = 0; i < Objects->size(); i++) { // TODO: re-add SEPERATE
        Region r(home); // TODO: one region per object or one for all objects? ("keep scope small")
        Object *o = Objects->collection[i];
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
              if (pMinStatus == ES_FAILED) {
                return ES_FAILED;
              }
            }
            F.resetRR(); // Reset RR position
            ExecStatus pMaxStatus = ES_FIX;
            if (!o->x[d].assigned()) {
              pMaxStatus = pruneMax(home, o, d, k, &F);
              if (pMaxStatus == ES_FAILED) {
                return ES_FAILED;
              }
            }
            if (pMinStatus == ES_NOFIX || pMaxStatus == ES_NOFIX) {
              nonfix = true; // We pruned a bound, not at fixpoint
            }
          }
        }

        if (!o->x.assigned()) {
          allfixed = false;
        }
      }
    }

    if (allfixed) {
      return home.ES_SUBSUMED(*this);
    }


    return ES_FIX;
  }

public:
  // Create propagator and initialize
  NonOverlapping(Home home, // Constructor for 2D
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
      o->l = ((Space&) home).alloc<int>(2);
      o->support_min = ((Space&) home).alloc<int*>(2);
      o->support_max = ((Space&) home).alloc<int*>(2);
      o->x = ViewArray<IntView>((Space&) home, 2);

      o->x[0] = x0[i];
      o->x[1] = y0[i];
      o->l[0] = w0[i];
      o->l[1] = h0[i];
      o->id = i;

      maxl[0] = std::max(maxl[0], o->l[0]);
      maxl[1] = std::max(maxl[0], o->l[1]);

      o->support_min[0] = ((Space&) home).alloc<int>(2);
      o->support_min[1] = ((Space&) home).alloc<int>(2);
      //o->support_min[1] = ((Space&) home).realloc<int>(o->support_min[0], 2, 2*2); // TODO: this doesn't seem worth it performance wise
      o->support_min[0][0] = o->x[0].min();
      o->support_min[0][1] = o->x[1].min();

      o->support_min[1][0] = o->x[0].min();
      o->support_min[1][1] = o->x[1].min();

      o->support_max[0] = ((Space&) home).alloc<int>(2);
      o->support_max[1] = ((Space&) home).alloc<int>(2);
      //o->support_max[1] = ((Space&) home).realloc<int>(o->support_max[0], 2, 2*2);
      o->support_max[0][0] = o->x[0].max();
      o->support_max[0][1] = o->x[1].max();

      o->support_max[1][0] = o->x[0].max();
      o->support_max[1][1] = o->x[1].max();

      o->rfre = ((Space&) home).alloc<int>(2);
      o->rfrb = ((Space&) home).alloc<int>(2);

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

    Int::IntView::schedule(home, *this, Int::ME_INT_DOM);
    home.notice(*this, AP_DISPOSE);
  }

  // Post no-overlap propagator
  static ExecStatus post(Home home,
                         ViewArray<IntView>& x, int w[],
                         ViewArray<IntView>& y, int h[]) {
    // Only if there is something to propagate
    if (x.size() > 1)
      (void) new (home) NonOverlapping(home,x,w,y,h);
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
  NonOverlapping(Space& home, bool share, NonOverlapping& p)
    : Propagator(home,share,p) {
    dimensions = p.dimensions;
    pivot = p.pivot;
    Objects = (OBJECTS*)((Space &) home).ralloc(sizeof(OBJECTS));
    new(Objects) OBJECTS((Space &) home);

    maxl = ((Space&) home).alloc<int>(2);
    
    for (int j = 0; j < dimensions; j++) {
      maxl[j] = p.maxl[j];
    }

    for (int i = 0; i < p.Objects->size(); i++) {
      Object *pObj = p.Objects->collection[i];
      Object *o = ((Space&) home).alloc<Object>(1);
      o->l = home.alloc<int>(dimensions);
      o->support_min = home.alloc<int*>(dimensions);
      o->support_max = home.alloc<int*>(dimensions);

      o->rfre = home.alloc<int>(dimensions);
      o->rfrb = home.alloc<int>(dimensions);

      o->d_size = home.alloc<int>(dimensions);

      o->x.update(home, share, pObj->x);
      for (int j = 0; j < dimensions; j++) {
        o->l[j] = pObj->l[j];
        o->support_min[j] = home.alloc<int>(dimensions);
        o->support_max[j] = home.alloc<int>(dimensions);
        for (int d = 0; d < dimensions; d++) {
          o->support_min[j][d] = pObj->support_min[j][d];
          o->support_max[j][d] = pObj->support_max[j][d];
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
    return new (home) NonOverlapping(home,share,*this);
  }
    
  // Return cost (defined as expensive quadratic)
  virtual PropCost cost(const Space&, const ModEventDelta&) const {
    return PropCost::quadratic(PropCost::HI,dimensions*Objects->size());
  }

  // TODO: why do we need this?
  virtual void reschedule(Space& home) {
    Int::IntView::schedule(home, *this, Int::ME_INT_DOM);
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
    return filter((Home) home, 2);
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
  if (NonOverlapping::post(home,vx,wc,vy,hc) != ES_OK)
    home.fail();
}


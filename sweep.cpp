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
  ForbiddenRegions(Region*);
};


ForbiddenRegions::ForbiddenRegions(Region *r) : collection(*r), length(0), RRpos(0) {}

FR* ForbiddenRegions::getRR() {
  return collection[(RRpos++ % length)];
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
  
  // Function for checking if o can possibly overlap f
  bool overlaps(FR *f, Object *o, int k) {
    for (int d = 0; d < k; d++) {
      if ((o->x[d].max() < f->min[d]) || (o->x[d].min() > f->max[d])) {
        return false;
      }
    }
    return true;
  }

  // Generates relative FRs to the object o and stores them in F
  void genOutBoxes(Home home, ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];

      if (!other->isSame(o)) { // For every object <> o
        FR *f = static_cast<Space&>(home).alloc<FR>(1);
        f->min = static_cast<Space&>(home).alloc<int>(dimensions);
        f->max = static_cast<Space&>(home).alloc<int>(dimensions);
        bool exists = true; // Assume f exists

        for (int d = 0; d < k; d++) { // For every dimension d
          const int min = other->x[d].max() - o->l[d] + 1;
          const int max = other->x[d].min() + other->l[d] - 1;
          if (min <= max) {
            f->min[d] = min;
            f->max[d] = max;
          } else {
            exists = false;
            break; // Break loop, other is not interfering with o
          }
        }

        if (!overlaps(f, o, k)) { // do not add if o does not possibly overlap with f
          exists = false;
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
  int mergeTest(FR *A, FR *B, int k) {
    int r = 0;
    for (int d = 0; d < k; d++) {
      if (A->min[d] <= B->min[d] && A->max[d] >= B->max[d]) {
        if (r == 2) {
          return 0;
        }
        r = 1;
      } else if (A->min[d] >= B->min[d] && A->max[d] <= B->max[d]) {
        if (r == 1) {
          return 0;
        }
        r = 2;
      } else {
        return 0;
      }
    }
    return r;
  }

  void genOutBoxesMerge(Space &home, ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    Region r(home);
    Support::StaticStack<FR*,Region> recent(r, 2); // merge stack for temporary storage of FRs (for merge checking)
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];

      if (!other->isSame(o)) { // For every other object <> o
        FR *f = (FR *) home.alloc<FR>(1);
        f->min = home.alloc<int>(dimensions);
        f->max = home.alloc<int>(dimensions);

        bool exists = true; // Assume f exists

        for (int d = 0; d < k; d++) { // For every dimension d
          const int min = other->x[d].max() - o->l[d] + 1;
          const int max = other->x[d].min() + other->l[d] - 1;
          if (min <= max) {
            f->min[d] = min;
            f->max[d] = max;
          } else {
            exists = false;
            break; // Break loop, other is not interfering with o
          }
        }

        if (!overlaps(f, o, k)) { // do not add if o does not overlap with f
          exists = false;
        }

        if (exists) {
          if (recent.entries() == 2) { // Check if previous 2 FRs can be merged
            FR *A = recent.pop();
            FR *B = recent.pop();
            int r = mergeTest(A, B, k);

            if (r == 0) { // no merge possible
              F->insert(B);
              recent.push(A);
            } else if (r == 1) { // B is inside A
              recent.push(A);
            } else { // A is inside B
              recent.push(B); 
            }
          }
          recent.push(f);
        }
      }
    }
    
    while (!recent.empty()) { // Make sure no FR is still on merge stack
      F->insert(recent.pop());
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
  ExecStatus pruneMin(Home home, ViewArray<IntView> o, int d, int k, ForbiddenRegions *F) {
    bool b = true;
    
    Region r(home); // Region for storage

    /* Init c and n*/
    int *c = r.alloc<int>(k); // sweep-point
    // init c
    for (int i = 0; i < k; i++) {
      c[i] = o[i].min();
    }

    int *n = r.alloc<int>(k); // jump vector
    // init vector
    for (int i = 0; i < k; i++) {
      n[i] = o[i].max() + 1;
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
        n[r] = o[r].max() + 1;    // Reset component of n after jumping
        if (c[r] <= o[r].max()) { // Jump target found?
          b = true;               // target found
          break;
        } else {
          c[r] = o[r].min();      // reset component of c, dimension r exhausted
        }
      }
      // Update currentF and check if it is infeasible still
      currentF = getFR(k, c, F);
    }

    if (b) { 
      /* If the sweep point arrived at a hole in the domain
       * adjust the sweep point in dimension d to next value in domain
       */
      // TODO: rerun pruning from that point since new point might be
      // infeasible? Due to foxpoint loop it might be redundant
      // however.
      if (!o[d].in(c[d])) { // TODO: can this be done more efficiently?
        for (IntVarValues i(o[d]); i(); ++i) {
          if (i.val() > c[d]) {
            c[d] = i.val();
            break;
          }
        }
      }

      ModEvent me = o[d].gq(home, c[d]); // prune o
      if (me_failed(me)) {
        return ES_FAILED;
      } else if (me_modified(me)) {
        return ES_NOFIX;
      } else {
        return ES_FIX;
      }
    } else {
      return ES_FAILED;
    }
  }

  // Function for pruning the upper bound of object o in dimension d
  ExecStatus pruneMax(Home home, ViewArray<IntView> o, int d, int k, ForbiddenRegions *F) {
    bool b = true;
    Region r(home);

    /* Init c and n*/
    int *c = r.alloc<int>(k); // sweep-point
    // init c
    for (int i = 0; i < k; i++) {
      c[i] = o[i].max();
    }    

    int *n = r.alloc<int>(k); // jump-vector
    // initvector
    for (int i = 0; i < k; i++) {
      n[i] = o[i].min() - 1;
    }

    // Get next FR
    FR *currentF = getFR(k, c, F);
    //bool infeasible = (currentF != NULL);

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
        n[r] = o[r].min() - 1;    // Reset component of n after jumping
        if (c[r] >= o[r].min()) { // Jump target found?
          b = true;               // target found
          break;
        } else {
          c[r] = o[r].max();      // reset component of c, dimension r exhausted
        }
      }

      // Update currentF and check if it is infeasible still
      currentF = getFR(k, c, F);
    }

    if (b) { 
      /* If the sweep point arrived at a hole in the domain
       * adjust the sweep point in dimension d to next value in domain
       */
      // TODO: rerun pruning from that point since new point might be
      // infeasible? Due to foxpoint loop it might be redundant
      // however.
      if (!o[d].in(c[d])) { 
        int prev = o[d].min() - 1;
        for (IntVarValues i(o[d]); i(); ++i) {
          if (i.val() >= c[d]) {
            c[d] = prev;
            break;
          }
          prev = i.val();
        }
      }

      ModEvent me = o[d].lq(home, c[d]); // prune o
      if (me_failed(me)) {
        return ES_FAILED;
      } else if (me_modified(me)) {
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
    Region r(home); // TODO: maybe new region for each object?

    bool nonfix = true;

    while (nonfix) {
      nonfix = false;
      for (int i = 0; i < pivot; i++) {
        Object *o = Objects->collection[i];
        
        ForbiddenRegions F(&r); 
        genOutBoxesMerge(home, &F, Objects, k, o);

        for (int d = 0; d < k; d++) {
          ExecStatus pMinStatus = pruneMin(home, o->x, d, k, &F);
          ExecStatus pMaxStatus = pruneMax(home, o->x, d, k, &F);
          if (pMinStatus == ES_FAILED || pMaxStatus == ES_FAILED) {
            return ES_FAILED;
          } else if (pMinStatus == ES_NOFIX || pMaxStatus == ES_NOFIX) {
            nonfix = true; // We pruned a bound, not at fixpoint
          }
        }
      }
    }

    // Check if any objects have become fixed
    for (int i = 0; i < pivot; i++) { // This pivot solution seems a bit slower compared to using fixed boolean attribute per object
        Object *o = Objects->collection[i];
        if (o->x.assigned()) {
          Objects->collection[i] = Objects->collection[pivot-1];
          Objects->collection[pivot-1] = o;
          --pivot; // Move pivot point to the left (we have one more fixed object)
          --i; // Need to check this position again, as it now contains other object
        }
    }

    if (pivot == 0) {
      //Objects->prettyPrint();
      return home.ES_SUBSUMED(*this);
    }


    return ES_FIX;
  }

public:
  // Create propagator and initialize
  NonOverlapping(Home home, // Constructor for 2D
                 ViewArray<IntView>& x0,int w0[],
                 ViewArray<IntView>& y0,int h0[])
    : Propagator(home)
  {
    Objects = (OBJECTS*)((Space &) home).ralloc(sizeof(OBJECTS));
    new(Objects) OBJECTS((Space &) home);

    // Create corresponding objects for the ViewArrays and arrays
    for (int i = 0; i < x0.size(); i++) {
      Object *o = ((Space&) home).alloc<Object>(1);
      o->l = ((Space&) home).alloc<int>(2);
      o->x = ViewArray<IntView>((Space&) home, 2);

      o->x[0] = x0[i];
      o->x[1] = y0[i];
      o->l[0] = w0[i];
      o->l[1] = h0[i];
      o->id = i;

      o->x.subscribe(home,*this,PC_INT_DOM);

      Objects->insert(o);
    }

    dimensions = 2;
    pivot = Objects->size(); // All objects left of pivot are non-fixed
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
    for (int i = 0; i < Objects->size(); i++) { 
      Objects->collection[i]->x.cancel(home, *this, PC_INT_DOM);
    }
    //home.ignore(* this,AP_DISPOSE);

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

    for (int i = 0; i < p.Objects->size(); i++) {
      Object *pObj = p.Objects->collection[i];
      Object *o = ((Space&) home).alloc<Object>(1);
      o->l = home.alloc<int>(dimensions);
      o->x.update(home, share, pObj->x);
      for (int j = 0; j < dimensions; j++) {
        o->l[j] = pObj->l[j];
      }
      o->id = pObj->id;
      Objects->insert(o);
    }

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
    for (int i = 0; i < Objects->size(); i++) {
      Objects->collection[i]->x.reschedule(home,*this,Int::PC_INT_DOM);
    }
  }
    
  // Perform propagation
  virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
    return filter((Home) home, 2);
    // ExecStatus fStatus = filter((Home) home, 2);// ? ES_FIX : ES_FAILED;

    // if (fStatus == ES_FAILED) {
    //   return ES_FAILED;
    // }

    // if (fStatus != ES_FIX) {
    //   return home.ES_SUBSUMED(*this);
    // } else {
    //   return fStatus;
    // }
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


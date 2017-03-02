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
#include <list> // TODO: remove

using namespace Gecode;
using namespace Gecode::Int;

class FR {
public:
  int *min; // origin per dimension
  int *max; // maximal per dimension
};

class Object {
public:
  ViewArray<IntView> x; // Object origins
  int *l; // Object lengths
  int id;
  bool fixed;
  
  bool isSame(Object *);

  Object() : fixed(false) {}
};

bool Object::isSame(Object *o) {
  return id == o->id;
}

class ForbiddenRegions {
private:
  int RRpos;
  int length;
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
public:
  Support::DynamicArray<Object*, Space> collection; 
  int size();
  int length;
  
  void insert(Object*);

  OBJECTS(Space& s) : collection(s), length(0) {};
};

int OBJECTS::size(void) {
  return length;
}

void OBJECTS::insert(Object *f) {
  collection[length] = f;
  ++length;
}

// The no-overlap propagator
class NonOverlapping : public Propagator {
protected:
  OBJECTS *Objects;
  //std::list<Object*> NonFixed;
  //std::list<Object*> Fixed;

  int dimensions;
  
  bool overlaps(FR *f, Object *o, int k) {
    for (int d = 0; d < k; d++) {
      if ((o->x[d].max() < f->min[d]) || (o->x[d].min() > f->max[d])) {
        return false;
      }
    }
    return true;
  }

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

  void genOutBoxesMerge(Home home, ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    Region r(home);
    Support::StaticStack<FR*,Region> recent(r, 2);
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];

      if (!other->isSame(o)) { // For every other object <> o
        FR *f = (FR *) static_cast<Space&>(home).alloc<FR>(1);//new FR(home, k);
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

        if (!overlaps(f, o, k)) { // do not add if o does not overlap with f
          exists = false;
        }

        if (exists) {
          if (recent.entries() == 2) {
            FR *A = recent.pop();
            FR *B = recent.pop();
            int r = mergeTest(A, B, k);

            if (r == 0) {
              F->insert(A);
              recent.push(B);
            } else if (r == 1) {
              recent.push(A);
            } else {
              recent.push(B);
            }
          }
          recent.push(f);
        }
      }
    }
    
    while (!recent.empty()) {
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

  FR *getFR(FR *f, int k, int *c, ForbiddenRegions *actrs) {
    for (int i = 0; i < actrs->size(); i++) {
      FR *actrsF = actrs->getRR();
      if (!isFeasible(actrsF, k, c)) { // forbidden region for c found
        return actrsF;
      }
    }
    return NULL;
  }

  virtual ExecStatus pruneMin(Home home, ViewArray<IntView> o, int d, int k, ForbiddenRegions *F) {
    bool b = true;
    ExecStatus pStatus = ES_FIX;
    
    Region r(home);

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

    FR *currentF = getFR(currentF, k, c, F);
    bool infeasible = (currentF != NULL);

    while (b and infeasible) {
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
      currentF = getFR(currentF, k, c, F);
      infeasible = (currentF != NULL);
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

      if (o[d].min() < c[d]) { // Is new min larger than old?
        pStatus = ES_NOFIX;
      } else {
        pStatus = ES_FIX;
      }

      GECODE_ME_CHECK(o[d].gq(home, c[d])); // prune minimum if target found
    } else {
      pStatus = ES_FAILED;
    }

    return pStatus;
  }

  virtual ExecStatus pruneMax(Home home, ViewArray<IntView> o, int d, int k, ForbiddenRegions *F) {
    bool b = true;
    ExecStatus pStatus = ES_FIX;

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

    FR *currentF = getFR(currentF, k, c, F);
    bool infeasible = (currentF != NULL);

    while (b and infeasible) {
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
      currentF = getFR(currentF, k, c, F);
      infeasible = (currentF != NULL);
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

      if (o[d].max() > c[d]) { // is new max smaller than old?
        pStatus = ES_NOFIX;
      } else {
        pStatus = ES_FIX; // We could not prune anything
      }

      GECODE_ME_CHECK(o[d].lq(home, c[d])); // prune minimum if target found
    } else {
      pStatus = ES_FAILED; // No feasible sweep point found, report failure
    }

    return pStatus;
  }

  // R is a collection of rectangles participating in the problem
  bool filter(Home home, int k) {
    Region r(home); // TODO: maybe new region for each object?
    //std::list<Object> *O = r.alloc< std::list<Object> >(1);
    //std::list<Object> *ToBeFiltered = static_cast<Space&>(home).alloc< std::list<Object> >(1);

    bool nonfix = true;

    //std::cout << "FILTER\n";
    while (nonfix) {
      nonfix = false;
      for (int i = 0; i < Objects->size(); i++) {
        Object *o = Objects->collection[i];

        if (o->fixed) { // TODO: this is likely a subpar SEPERATE-implementation
          continue;          
        }
        
        ForbiddenRegions F(&r); 
        genOutBoxesMerge(home, &F, Objects, k, o);

        for (int d = 0; d < k; d++) {
          ExecStatus pMinStatus = pruneMin(home, o->x, d, k, &F);
          ExecStatus pMaxStatus = pruneMax(home, o->x, d, k, &F);
          if (pMinStatus == ES_FAILED || pMaxStatus == ES_FAILED) {
            return false;
          } else if (pMinStatus == ES_NOFIX || pMaxStatus == ES_NOFIX) {
            nonfix = true; // We pruned a bound, not at fixpoint
          }
        }

      }
    }

    // Check if any objects have become fixed
    for (int i = 0; i < Objects->size(); i++) {
        Object *o = Objects->collection[i];
        if (o->x.assigned() && !o->fixed) {
          o->x.cancel(home, *this, PC_INT_DOM);
          o->fixed = true;
        }
    }

    return true;
  }

public:
  // Create propagator and initialize
  NonOverlapping(Home home,
                 ViewArray<IntView>& x0,int w0[],
                 ViewArray<IntView>& y0,int h0[])
    : Propagator(home)
  {
    Objects = (OBJECTS*)((Space &) home).ralloc(sizeof(OBJECTS));
    new(Objects) OBJECTS((Space &) home);

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

  // Dispose propagator and return its sizestd::list<Object*>
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
      o->fixed = pObj->fixed;
      Objects->insert(o);
    }

  }
  // Create copy during cloning
  virtual Propagator* copy(Space& home, bool share) {
    return new (home) NonOverlapping(home,share,*this);
  }
    
  // Return cost (defined as cheap quadratic)
  virtual PropCost cost(const Space&, const ModEventDelta&) const {
    return PropCost::quadratic(PropCost::LO,2*((int) 2));
  }

  // TODO: why do we need this?
  virtual void reschedule(Space& home) {
    for (int i = 0; i < Objects->size(); i++) {
      Objects->collection[i]->x.reschedule(home,*this,Int::PC_INT_DOM);
    }
  }
    
  // Perform propagation
  virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
    //filter((Home) home, 2, Objects);       
        
    ExecStatus fStatus = filter((Home) home, 2) ? ES_FIX : ES_FAILED;

    if (fStatus == ES_FAILED) {
      return ES_FAILED;
    }

    bool subsumed = true;
    int i,j;

    // TODO: this makes things a tad slower, figure out how to perform subsumption check on the go?
    for (int f = 0; f < Objects->size(); f++) {
      Object *first = Objects->collection[f];
      for (int s = 0; s < Objects->size(); s++) {
        Object *second = Objects->collection[s];
        if (first->id == second->id) { // Don't run on a pair of identical rectangles
          continue;
        }

        // // x overlap
        // if (first->x[0].min()+first->l[0] > second->x[0].max() && first->x[0].max() < second->x[0].min()+second->l[0]) {
        //   // i cannot be above j
        //   if (first->x[1].max() < second->x[1].min()+second->l[1]) {
        //     GECODE_ME_CHECK(first->x[1].lq(home, second->x[1].max()-first->l[1]));
        //     GECODE_ME_CHECK(second->x[1].gq(home, first->x[1].min()+first->l[1]));
        //   }
        // }

        // //y-overlap
        // if (first->x[1].min()+first->l[1] > second->x[1].max() && first->x[1].max() < second->x[1].min()+second->l[1]) {
        //   // i cannot be to the right of j                    
        //   if (first->x[0].max() < second->x[0].min()+second->l[0]) {
        //     GECODE_ME_CHECK(first->x[0].lq(home, second->x[0].max()-first->l[0]));
        //     GECODE_ME_CHECK(second->x[0].gq(home, first->x[0].min()+first->l[0]));
        //   } 
        // }

          //Subsumption check
          //Potential overlap on x-axis
        if (first->x[0].max() + first->l[0] > second->x[0].min() && second->x[0].min() >= first->x[0].min()) {
          // Potential overlap on y-axis
          if(first->x[1].max() + first->l[1] > second->x[1].min() && second->x[1].min() >= first->x[1].min()) {
            subsumed = false; // Potential overlap on both axis, can't subsume
          }
        }
      }
    }

    if (subsumed) {
      //std::cout << "Subsumed\n";
      return home.ES_SUBSUMED(*this);
    } else {
      return fStatus;
      //return ES_NOFIX;
    }
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

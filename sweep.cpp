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

/* The Dim struct represents a minimum and maximum value in a given
   dimension */
struct Dim {
  int min;
  int max;
};

/* The FR struct represents forbidden regions consisting of k
   dimensions, each dimension containing min and max of the FR */
struct FR {
  Dim dim[];
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

  // Reference counter used in ENABLE optimisation, if > 0, then the
  // object can be skipped during relative FR generation
  int skippable;
  
  bool isSame(Object *);
};

forceinline
bool Object::isSame(Object *o) {
  return id == o->id;
}

/* The ForbiddenRegions class is used for storing forbidden regions
   and retrieving them using RR scheduling */
class ForbiddenRegions {
private:
  int RRpos;          // Current RR position in GetFR
  int length;         // Amount of FRs in collection
  int dimensions;     // Number of dimensions used
public:
  FR* collection;     // Collection stores FRs in place
public:
  int size();         // Gives the amount of FRs in collection
  void insert(FR*);   // Insert FR at top of stack (by copying contents)
  FR* getRR();        // Get FR based on RR scheduling
  void resetRR();     // Reset RR scheduling
  void prettyPrint(); // Prettyprinting the FRs
  FR* get(int);       // Get FR at given index
  FR* getTop();       // Get FR on top of FR stack
  void incTop();      // Increment the top of the stack
  void decTop();      // Decrement the top of the stack
  ForbiddenRegions(Region &r, int, int);
};


ForbiddenRegions::ForbiddenRegions(Region &r, int k, int n) : length(0), RRpos(0), dimensions(k) {
  collection = (FR *) r.ralloc((sizeof(FR) + sizeof(Dim)*k)*n);
}

forceinline FR *ForbiddenRegions::getTop() {
  return (FR*) ((char *)collection + (sizeof(FR) + sizeof(Dim)*dimensions)*length);
}

forceinline FR *ForbiddenRegions::get(int idx) {
  return (FR*) ((char *)collection + (sizeof(FR) + sizeof(Dim)*dimensions)*idx);
}

forceinline void ForbiddenRegions::incTop() {
  length++;
}

forceinline void ForbiddenRegions::decTop() {
  length--;
}

forceinline void ForbiddenRegions::prettyPrint() {
  std::cout << "----------ForbiddenRegions----------\n";
  for (int i = 0; i < size(); i++) {
    FR *f = get(i);
    std::cout << "(";
    
    for (int d = 0; d < dimensions; d++) {
      std::cout << f->dim[d].min << ".." << f->dim[d].max;

      if (not d + 1 == dimensions) {
        std::cout << ",";
      }
    }
    std::cout << ")\n";
  }
  std::cout << "------------------------------------\n";
}

forceinline FR* ForbiddenRegions::getRR() {
  return get(RRpos++ % length);
}

forceinline void ForbiddenRegions::resetRR() {
  RRpos = 0;
}

forceinline int ForbiddenRegions::size(void) {
  return length;
}

forceinline void ForbiddenRegions::insert(FR *f) {
  for (int i = 0; i < dimensions; i++) {
    get(length)->dim[i].min = f->dim[i].min;
    get(length)->dim[i].max = f->dim[i].max;
  }
  incTop();
}

/* The OBJECTS class, used for storing Objects */
class OBJECTS {
  int length;
public:
  Object** collection; // Collection stores pointers to Objects
  int size();
  void insert(Object*, int i);
  void prettyPrint();
  void removeLast();

  OBJECTS(Space &h, int l);
};

OBJECTS::OBJECTS(Space &h, int l) {
  length = l;
  collection = h.alloc<Object*>(length);
}

forceinline void OBJECTS::prettyPrint() {
  for (int i = 0; i < size(); i++) {
    Object *o = collection[i];
    std::cout << o->x << "\n";
  }
}

forceinline int OBJECTS::size(void) {
  return length;
}

forceinline void OBJECTS::removeLast(void) {
  length--;
}

forceinline void OBJECTS::insert(Object *f, int i) {
  collection[i] = f;
}

// The diffn propagator
class Diffn : public Propagator {
protected:
  OBJECTS *Objects; // Objects being filtered
  int dimensions; // Number of dimensions of the problem

  FR *B; // Bounding box for INCREMENTAL optimsation
  bool inFilter; // Flag indacting if a bound was changed by filter or not

  int pointerSize;
  int *Pointers; // Pointers[id] denotes the index in Objects for object id
  int pivot;

  int *maxl; // The recorded maximum length of all objects

  class ViewAdvisor : public Advisor {
  public:
    //Int::IntView x;
    int dim; // What dimension in the corresponding viewarray is this advisor responsible for?
    int i; // The corresponding object
    bool skippable; // ENABLE optimisation: is the corresponding object skippable in this dimension?
    ViewAdvisor(Space& home, Propagator& p, 
                Council<ViewAdvisor>& c, int d, int idx) 
      : Advisor(home,p,c), dim(d), i(idx), skippable(true) {
    }
    ViewAdvisor(Space& home, Propagator& p, 
                Council<ViewAdvisor>& c, int d, int idx, bool skipped) 
      : Advisor(home,p,c), dim(d), i(idx), skippable(skipped) {
    }
    ViewAdvisor(Space& home, bool share, ViewAdvisor& a)
      : Advisor(home,share,a), dim(a.dim), i(a.i), skippable(a.skippable) {
    }
    void dispose(Space& home, Council<ViewAdvisor>& c) {
      Advisor::dispose(home,c);
    }
  };
  
  Council<ViewAdvisor> c;
  
  // True iff (min..max) and o overlap in dimension d
  bool overlaps(int min, int max, Object *o, int d) {
    if ((o->x[d].max() < min) || (o->x[d].min() > max)) {
      return false;
    }
    return true;
  }

  // True iff o cannot intersect the boundingBox
  bool cantOverlap(FR *boundingBox, Object *o, int k) {
    for (int i = 0; i < k; i++) {
      if ((o->x[i].max() + o->l[i] - 1 < boundingBox->dim[i].min) || (o->x[i].min() > boundingBox->dim[i].max)) {
        return true;
      }
    }
    return false;
  }

  // Generates relative FRs to the object o and stores them in F
  forceinline void genOutBoxes(ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    int top = F->size();
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];

      if (!other->isSame(o)) { // For every object <> o
        if (other->skippable > 0) { 
          continue;
        }

        FR *f = F->get(top);
        bool exists = true; // Assume f exists

        for (int d = 0; d < k; d++) { // For every dimension d
          const int min = other->rfrb[d] - o->l[d];
          const int max = other->rfre[d];
          if ((min <= max) && overlaps(min, max, o, d)) {
            f->dim[d].min = min;
            f->dim[d].max = max;
          } else {
            exists = false;
            break; // Break loop, other is not interfering with o
          }
        }

        // Only add forbidden region if it exists
        if (exists) {
          F->incTop();
          top++;
        }
      }
    }
  }
  
    /* Tries to coalesce f0 and f1.
   * Return cases:
   * 0: coalescing not possible
   * 1: f0 includes f1 (completely overlaps) or one of the two can be extended in one dimension to represent both FRs,
   *    f0 is then changed in place to represent both.
   * 2: f1 includes f2 (completely overlaps)
   */
  forceinline
  int coalesce(FR *f0, FR *f1, int k) {
    // TODO: introduce enumerator
    int trend = 0; // What case did the previous dimensions fullfil?
    int e = 0; // In case 3, what dimension can be extended?

    for (int d = 0; d < k; d++) {
      if (f0->dim[d].max + 1 < f1->dim[d].min || f0->dim[d].min > f1->dim[d].max + 1) { /* Do not coalesce */
        return 0;
      } else if (f0->dim[d].min == f1->dim[d].min && f0->dim[d].max == f1->dim[d].max) { /* Equal in dimension d */
        // no-op
      } else if (f0->dim[d].min <= f1->dim[d].min && f0->dim[d].max >= f1->dim[d].max) { /* f1 \subset f0*/
        if (trend == 0 || trend == 1) {
          trend = 1;
        } else {
          return 0; // The FRs cannot be coalesced
        }
      } else if (f0->dim[d].min >= f1->dim[d].min && f0->dim[d].max <= f1->dim[d].max) { /* f0 \subset f1*/
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
      f0->dim[e].min = std::min(f0->dim[e].min, f1->dim[e].min);
      f0->dim[e].max = std::max(f0->dim[e].max, f1->dim[e].max);
      return 1;
    }
  }

  /* Generates forbidden regions for object o given objects in O, merges forbidden regions where possible */
  forceinline void genOutBoxesMerge(ForbiddenRegions *F, OBJECTS *O, int k, Object *o) {
    for (int i = 0; i < O->size(); i++) {
      Object *other = O->collection[i];
      
      if (!other->isSame(o)) { // For every other object <> o
        if (other->skippable > 0) { 
          continue;
        }

        FR *f = F->getTop();

        bool exists = true; // Assume f exists

        for (int d = 0; d < k; d++) { // For every dimension d
          const int min = other->rfrb[d] - o->l[d];
          const int max = other->rfre[d];
          if ((min <= max) && overlaps(min, max, o, d)) {
            f->dim[d].min = min;
            f->dim[d].max = max;
          } else {
            exists = false;
            break; // Break loop, other is not interfering with o
          }
        }

        if (exists) {
          while (exists && F->size() > 0) {
            FR *f0 = F->get(F->size()-1); // Get the FR one step from the top
            switch (coalesce(f0, f, k)) { /* Try to coalesce the FRs */
            case 0:  /* No coalescing possible */
              // no-op, f0 stays in place
              exists = false;
              break;
            case 1: /* f0 subsumes f */
              F->decTop();
              f = f0;
              break;
            case 2: /* f subsumes f0 */
              for (int d = 0; d < k; d++) {
                f0->dim[d].min = f->dim[d].min;
                f0->dim[d].max = f->dim[d].max;
              }
              f = f0;
              F->decTop();
              break;
            }
          }
          F->incTop();
        }
      }
    }
  }

  // Returns true iff c is feasible according to outbox
  forceinline bool isFeasible(FR *outbox, int k, int *c) {
    for (int d = 0; d < k; d++) {
      if ((c[d] < outbox->dim[d].min) || (c[d] > outbox->dim[d].max)) {
        return true;
      }
    }
    return false;
  }

  // Returns next infeasible FR if one exists, NULL otherwise
  forceinline FR *getFR(int k, int *c, ForbiddenRegions *actrs) {
    for (int i = 0; i < actrs->size(); i++) {
      FR *f = actrs->getRR(); // Using RR for getting the next FR
      if (!isFeasible(f, k, c)) { // forbidden region for c found
        return f; // f is an infeasible FR
      }
    }
    return NULL; // return NULL if no infeasible FR was found
  }

  // Adjust the sweep-point based on jump vector information
  forceinline bool adjustUp(Object *o, int *c, int *n, int d, int k) {
    bool b = false;// No new point to jump to yet (assume failure)

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

    return b;
  }

  // Function for pruning the lower bound of object o in dimension d
  forceinline ExecStatus pruneMin(Home home, Object *o, int d, int k, ForbiddenRegions *F) {
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
      for (int i = 0; i < k; i++) {
        n[i] = std::min(n[i], currentF->dim[i].max + 1);
      }
      
      b = adjustUp(o, c, n, d, k);

      // Update currentF and check if it is infeasible still
      currentF = getFR(k, c, F);
    }

    if (b) { 
      /* SUPPORT optimisation */
      for (int j = 0; j < k; j++) {
        o->support_min[d*k+j] = c[j];
      }

      ModEvent me = o->x[d].gq(home, c[d]); // prune o
      if (me_modified(me)) {
        return ES_NOFIX;
      } else {
        return ES_FIX;
      }
    } else {
      return ES_FAILED;
    }
  }

  forceinline bool adjustDown(Object *o, int *c, int *n, int d, int k) {
    bool b = false; // No new point to jump to yet (assume failure)

    for (int j = k - 1; j >= 0; j--) {
      int r = (j + d) % k;
      c[r] = n[r];                 // Use n vector to jump
      n[r] = o->x[r].min() - 1;    // Reset component of n after jumping
      if (c[r] >= o->x[r].min()) { // Jump target found?
        b = true;                  // target found
        break;
      } else {
        c[r] = o->x[r].max();      // reset component of c, dimension r exhausted
      }
    }

    return b;
  }


  // Function for pruning the upper bound of object o in dimension d
  forceinline ExecStatus pruneMax(Home home, Object *o, int d, int k, ForbiddenRegions *F) {
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
      for (int i = 0; i < k; i++) {
        n[i] = std::max(n[i], currentF->dim[i].min - 1);
      }

      b = adjustDown(o, c, n, d, k);

      // Update currentF and check if it is infeasible still
      currentF = getFR(k, c, F);
    }

    if (b) { 
      /* SUPPORT optimisation */
      for (int j = 0; j < k; j++) {
        o->support_max[d*k+j] = c[j];
      }

      ModEvent me = o->x[d].lq(home, c[d]); // prune o
      if (me_modified(me)) {
        return ES_NOFIX;
      } else {
        return ES_FIX;
      }
    } else {
      return ES_FAILED; // No feasible sweep point found, report failure
    } 
  }

  // R is a collection of rectangles participating in the problem
  forceinline ExecStatus filter(Space &home, int k) {
    inFilter = true; // Entering filter - future bound changes are made by filter
    bool nonfix = true;

    // Bounding box for temporary storage of B
    FR *internalB = (FR *) home.ralloc(sizeof(FR) + sizeof(Dim)*k);

    while (nonfix) {
      nonfix = false;

      for (int j = 0; j < k; j++) {
        // Move values from B to internalB (so that B can be populated by internal events)
        internalB->dim[j].min = B->dim[j].min;
        internalB->dim[j].max = B->dim[j].max;
        
        // Reset bounding box B to consider internal events within fixpoint loop
        B->dim[j].min = Gecode::Int::Limits::infinity;
        B->dim[j].max = Gecode::Int::Limits::min;
      }

      for (int i = 0; i <= pivot; i++) {
        Object *o = Objects->collection[i];

        if (cantOverlap(internalB, o, k)) { // Consider external events
          continue; 
        }

        Region r(home);
        ForbiddenRegions F(r, k, Objects->size()-1); 
        genOutBoxes(&F, Objects, k, o);
        
        if (o->x.assigned()) {
          if (F.size() > 0) { // If a FR exists, then o->x must be infeasible
            return ES_FAILED;
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

        if (o->x.assigned()) {
          // swap as o has become fixed (and checked)
          Object *pivotObject = Objects->collection[pivot];
          Objects->collection[i] = pivotObject;
          Objects->collection[pivot] = o;
          Pointers[pivotObject->id] = i; // Advisors for pivotObject now search at position i for pivotObject
          pivot--; // Decrement pivot as one more object is now fixed
          i--; // Decrement so we don't miss swapped object
        } 
      }
    }

    // Free space for internalB so that it may be reused
    home.rfree(internalB, sizeof(FR) + sizeof(Dim)*k);

    Region r(home);
    // activeBox is a bounding box for all non-fixed objects
    FR *activeBox = (FR *) r.ralloc((sizeof(FR) + sizeof(Dim)*k));
    
    for (int j = 0; j < k; j++) {
      activeBox->dim[j].min = Gecode::Int::Limits::infinity;
      activeBox->dim[j].max = Gecode::Int::Limits::min;
    }

    for (int i = 0; i <= pivot; i++) {
      Object *o = Objects->collection[i];

      for (int j = 0; j < k; j++) {
        activeBox->dim[j].min = std::min(activeBox->dim[j].min, o->x[j].min());
        activeBox->dim[j].max = std::max(activeBox->dim[j].max, o->x[j].max() + o->l[j] - 1);
      }
    }

    for (int i = pivot + 1; i < Objects->size(); i++) {
      Object *o = Objects->collection[i];

      if (cantOverlap(activeBox, o, k)) {
        // Remove o from Objects
        Objects->collection[i] = Objects->collection[Objects->size() - 1];
        Objects->removeLast();
        --i; // Must check new i!
      }
    }

    // Reset bounding box
    for (int j = 0; j < k; j++) {
      B->dim[j].min = Gecode::Int::Limits::infinity;
      B->dim[j].max = Gecode::Int::Limits::min;
    }

    inFilter = false; // Exiting filter - future bound changes are made externally

    // If all objects are fixed and we have not failed, we can subsume
    if (pivot < 0) {
      return home.ES_SUBSUMED(*this);
    }

    return ES_FIX; // Filter is a fixpoint loop, thus report at fixpoint
  }

public:
  // Create propagator and initialize
  Diffn(Home home, // Constructor for 2D
        const IntVarArgs& x0,int w0[],
        const IntVarArgs& y0,int h0[])
    : Propagator(home), c(home)
  {
    Objects = (OBJECTS*)((Space &) home).ralloc(sizeof(OBJECTS));
    new(Objects) OBJECTS((Space &) home, x0.size());

    inFilter = false;

    pivot = x0.size()-1; // Everything to the right of pivot is fixed
    pointerSize = x0.size();
    Pointers = (int*)((Space &) home).ralloc(sizeof(int)*pointerSize);

    maxl = ((Space&) home).alloc<int>(2);
    maxl[0] = -1;
    maxl[1] = -1;

    // Create corresponding objects for the ViewArrays and arrays
    for (int i = 0; i < x0.size(); i++) {
      Object *o = ((Space&) home).alloc<Object>(1);

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
      maxl[1] = std::max(maxl[1], o->l[1]);

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

      o->skippable = 2; // Assume skippable in both dimensions

      Objects->insert(o, i);
      Pointers[i] = i;
    }

    /* Need to calculate skippable here since maxl is not calculated fully above */
    for (int i = 0; i < x0.size(); i++) {
      Object *o = Objects->collection[i];
      bool tmp[2] = {true, true}; // tmp[i] = true iff o is skippable in dimension i

      for (int j = 0; j < 2; j++) {
        if (not (o->rfrb[j] - o->rfre[j] > maxl[j])) {
          o->skippable--;
          tmp[j] = false;
        }
      }

      o->x[0].subscribe(home,*new (home) ViewAdvisor(home,*this,c,0,i,tmp[0]));
      o->x[1].subscribe(home,*new (home) ViewAdvisor(home,*this,c,1,i,tmp[1]));
    }

    dimensions = 2;

    B = (FR *) ((Space &) home).ralloc((sizeof(FR) + sizeof(Dim)*dimensions));

    for (int i = 0; i < dimensions; i++) {
      // Make sure every rectangle is checked at first filtering
      B->dim[i].min = Gecode::Int::Limits::min;
      B->dim[i].max = Gecode::Int::Limits::infinity;
    }

    Int::IntView::schedule(home, *this, Int::ME_INT_BND); // Schedule the propagator
    home.notice(*this, AP_DISPOSE); // Make sure dispose function is called on Space destruction
  }

  // Post no-overlap propagator
  static ExecStatus post(Home home,
                         const IntVarArgs& x, int w[],
                         const IntVarArgs& y, int h[]) {
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
    new(Objects) OBJECTS((Space &) home, p.Objects->size());

    inFilter = false;

    pivot = p.pivot;
    pointerSize = p.pointerSize;
    Pointers = (int*)((Space &) home).ralloc(sizeof(int)*p.pointerSize);

    B = (FR *) home.ralloc((sizeof(FR) + sizeof(Dim)*dimensions));

    for (int j = 0; j < dimensions; j++) {
      // Reset bounding box
      B->dim[j].min = Gecode::Int::Limits::infinity;
      B->dim[j].max = Gecode::Int::Limits::min;
    }

    maxl = ((Space&) home).alloc<int>(dimensions);
    
    for (int j = 0; j < dimensions; j++) {
      maxl[j] = p.maxl[j];
    }

    for (int i = 0; i < p.pointerSize; i++) {
      Pointers[i] = p.Pointers[i];
    }


    for (int i = 0; i < p.Objects->size(); i++) {
      Object *pObj = p.Objects->collection[i];
      Object *o = ((Space&) home).alloc<Object>(1);

      o->l = home.alloc<int>(dimensions*3); // Make sure memory block fits 2 more arrays of identical size
      o->rfre = &(o->l[dimensions]);
      o->rfrb = &(o->rfre[dimensions]);

      for (int j = 0; j < dimensions; j++) {
        o->l[j] = pObj->l[j];
        o->rfre[j] = pObj->rfre[j];
        o->rfrb[j] = pObj->rfrb[j];
      }

      o->x.update(home, share, pObj->x);
      o->skippable = pObj->skippable;
      o->id = pObj->id;
      Objects->insert(o, i);
    }

    // Only copy support if the object is not fixed
    for (int i = 0; i <= pivot; i++) {
      Object *pObj = p.Objects->collection[i];
      Object *o = Objects->collection[i];
      
      o->support_min = (int *) home.ralloc(sizeof(int) * dimensions * dimensions); // 1D representation of matrix
      o->support_max = (int *) home.ralloc(sizeof(int) * dimensions * dimensions); // 1D representation of matrix
      
      for (int j = 0; j < dimensions; j++) {
        for (int d = 0; d < dimensions; d++) {
          o->support_min[j * dimensions + d] = pObj->support_min[j * dimensions + d];
          o->support_max[j * dimensions + d] = pObj->support_max[j * dimensions + d];
        }
      }
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

  // Advise function, scheduled whenever its corresponding view changes
  virtual ExecStatus advise(Space& home, Advisor& _a, const Delta& d) {
    ModEvent me = IntView::modevent(d);

    /* Only update values if a bound was changed */
    if (me == ME_INT_BND || me == ME_INT_VAL) {
      ViewAdvisor& a = static_cast<ViewAdvisor&>(_a);
      int dim = a.dim;
      Object *o = Objects->collection[Pointers[a.i]];
      
      // update values since view changed
      o->rfrb[dim] = o->x[dim].max() + 1; 
      o->rfre[dim] = o->x[dim].min() + o->l[dim] - 1;
    
      if (a.skippable) {
        if (not (o->rfrb[dim] - o->rfre[dim] > maxl[dim])) {
          o->skippable--;
          a.skippable = false; // Don't decrement skippable counter more than once per dimension
        }
      }

      // Update bounding box B
      for (int j = 0; j < dimensions; j++) {
        B->dim[j].min = std::min(B->dim[j].min, o->x[j].min());
        B->dim[j].max = std::max(B->dim[j].max, o->x[j].max() + o->l[j] - 1);
      }
      
      if (me == ME_INT_VAL) {
        // Dispose advisor as its view became fixed
        return (inFilter == false) ? home.ES_NOFIX_DISPOSE(c, a) : home.ES_FIX_DISPOSE(c, a);
      } else {
        // Only schedule if view was modified externally
        return (inFilter == false) ? ES_NOFIX : ES_FIX;
      }
    } else {
      return ES_FIX; // Don't schedule if no bound was changed
    }
  }
    
  // Perform propagation
  virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
    return filter(home, dimensions); // TODO: fixing dimensions = 2 here gives performance boost
  }
 
};

/*
 * Post the constraint that the rectangles defined by the coordinates
 * x and y and width w and height h do not overlap.
 *
 * This is the function that you will call from your model. The best
 * is to paste the entire file into your model.
 */
void diffn(Home home,
                    const IntVarArgs& x, const IntArgs& w,
                    const IntVarArgs& y, const IntArgs& h) {
  // Check whether the arguments make sense
  if ((x.size() != y.size()) || (x.size() != w.size()) ||
      (y.size() != h.size()))
    throw ArgumentSizeMismatch("nooverlap");
  // Never post a propagator in a failed space
  if (home.failed()) return;
  // Set up arrays (allocated in home) for width and height and initialize
  int* wc = static_cast<Space&>(home).alloc<int>(x.size());
  int* hc = static_cast<Space&>(home).alloc<int>(y.size());
  for (int i=x.size(); i--; ) {
    wc[i]=w[i]; hc[i]=h[i];
  }
  // If posting failed, fail space
  if (Diffn::post(home,x,wc,y,hc) != ES_OK)
    home.fail();
}


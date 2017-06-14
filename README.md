# Sweep-Based Propagator for Diffn

An implementation of a sweep-based propagator for the [Diffn](http://web.emn.fr/x-info/sdemasse/gccat/Cdiffn.html) constraint in the copying constraint programming solver [Gecode](http://www.gecode.org/). Builds on the ideas of Beldiceanu and Carlsson from their paper [Sweep as a generic pruning technique applied to the non-overlapping rectangles constraint (2001)](https://link.springer.com/chapter/10.1007%2F3-540-45578-7_26).

# Usage

## 1D

Include the file `sweep.cpp` in a constraint programming model written in Gecode and call `diffn(home, x, w)`, where:

* `x` : IntVarArgs represent the origin of the lines (defined as the lexicographically smallest point)
* `w` : IntArgs represent the size of the lines


## 2D

Include the file `sweep.cpp` in a constraint programming model written in Gecode and call `diffn(home, x, w, y, h)`, where:

* `x` : IntVarArgs represent the origin in the *x*-axis of the rectangles (defined as the lexicographically smallest point)
* `w` : IntArgs represent the widths of the rectangles
* `y` : IntVarArgs represent the origin in the *y*-axis of the rectangles (defined as the lexicographically smallest point)
* `h` : IntArgs represent the heights of the rectangles

## 3D


Include the file `sweep.cpp` in a constraint programming model written in Gecode and call `diffn(home, x, w, y, h, z, l)`, where:

* `x` : IntVarArgs represent the origin in the *x*-axis of the cuboids (defined as the lexicographically smallest point)
* `w` : IntArgs represent the widths of the rectangles
* `y` : IntVarArgs represent the origin in the *y*-axis of the cuboids (defined as the lexicographically smallest point)
* `h` : IntArgs represent the heights of the rectangles
* `z` : IntVarArgs represent the origin in the *z*-axis of the cuboids (defined as the lexicographically smallest point)
* `l` : IntArgs represent the length in the *z*-dimension of the rectangles

## *k*-D

Include the file `sweep.cpp` in a constraint programming model written in Gecode and call `diffn(home, x, l, k)`, where:

* `x`: IntVarArgs represent the origin coordinates of the hyperrectangle (defined as the lexicographically smallest point). Here, hyperrectangle object *i* is represented by the coordinates `x[i*k]`, `x[i*k+1]`, .., `x[i*k+k]`. Where `x[i*k]` represents *i*'s coordinate in dimension 0 (or x), `x[i*k+1]` its coordinate in dimension 1 (or y) etc.
* `l`: IntArgs represent the side-length of the hyperrectangles. Here, the side-lengths of a hyperrectangle object *i* is represented by the values `l[i*k]`, `l[i*k+1]`, .., `i[i*k+k]`. 
* `k`: Int represents the dimensionality of the hyperrectangles.

# Report

A MSc thesis describing and evaluating the propagator is currently awaiting publication. This entry in the readme will later point you to it.

# Notes

The propagator currently does not support defining the hyperrectnagles' side-lengths as variables.

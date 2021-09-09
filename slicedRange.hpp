// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#ifndef SLICEDRANGE_HPP
#define SLICEDRANGE_HPP

#include <cmath>

template <typename T> class slicedRange {
public:
  slicedRange() {
    vmin = 0;
    n = 0;
    scale = 1.0;
  }

  slicedRange(T Vmin, T Vmax, int N) : vmin(Vmin), n(N) { scale = (T)n / (Vmax - Vmin); }

  void set(T Vmin, T Vmax, int N) {
    vmin = Vmin;
    n = N;
    scale = (T)n / (Vmax - Vmin);
  }

  void set_MinMaxNb(T Vmin, T Vmax, int N) { set(Vmin, Vmax, N); }
  void set_MinWidthNb(T Vmin, T W, int N) { set(Vmin, Vmin + (double)N * W, N); }

  // Returns a negative value if not in range
  int getID(T value) {
    int i = (int)floor(scale * (value - vmin));
    if (i >= n)
      i = -i;
    return i;
  }

  T getStep() const { return (1.0f / scale); }
  int getNumberOfSlices() const { return n; }
  T getLeftValue() const { return vmin; }
  T getRightValue() const { return vmin + (double)n / scale; }

protected:
  T vmin;  // left value
  int n;   // Number of slices
  T scale; // inverse of the bin width
};

#endif /* end of include guard: SLICEDRANGE_HPP */

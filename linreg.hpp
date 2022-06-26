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

#ifndef LINREG_HPP
#define LINREG_HPP

#include <cmath>
#include <iostream>
#include <vector>

class linreg {
public:
  static linreg *get() {
    static linreg instance;
    return &instance;
  }

  double orig;
  double slope;
  double err;

  /// @brief Make a linear regression of a set of (x, y) pairs
  ///        it computes orig, slope and err so that the line orig x + +slope is the linear regression of data x and y.
  ///        err is the coefficient of correlation.
  void run(std::vector<double> &x, std::vector<double> &y) {
    if (x.size() != y.size()) {
      std::cerr << "@linreg, x and y must have the same size" << std::endl;
      return;
    }
    int N = static_cast<int>(x.size());
    if (N == 0)
      return;
    double invN = 1.0f / (double)N;

    double xmean = 0.0;
    double ymean = 0.0;
    for (int i = 0; i < N; i++) {
      xmean += x[i];
      ymean += y[i];
    }
    xmean *= invN;
    ymean *= invN;

    double Sx = 0.0;
    double Sy = 0.0;
    double Sxy = 0.0;
    for (int i = 0; i < N; i++) {
      double dx = x[i] - xmean;
      double dy = y[i] - ymean;
      Sx += dx * dx;
      Sy += dy * dy;
      Sxy += dx * dy;
    }

    orig = Sxy / Sx; // suppose that Sx cannot be zero
    slope = ymean - orig * xmean;
    if (Sy == 0.0) {
      err = 1.0;
    } else {
      err = Sxy / sqrt(Sx * Sy);
    }
  }

private:
  linreg() {}
};

#endif /* end of include guard: LINREG_HPP */

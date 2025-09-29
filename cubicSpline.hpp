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

#ifndef CUBICSPLINE_HPP
#define CUBICSPLINE_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

struct SplineSet {
  double a;
  double b;
  double c;
  double d;
  double x;
};

/**
 * Computes the coefficients of a cubic spline interpolation for a given set of
 * data points.
 *
 * @param x A vector of doubles representing the x-coordinates of the data points.
 *          The vector must be sorted in ascending order.
 * @param y A vector of doubles representing the y-coordinates of the data points.
 *          Must be the same size as x.
 * @return A vector of SplineSet objects, each containing the coefficients
 *         'a', 'b', 'c', 'd', and 'x' for a segment of the spline.
 *
 * This function calculates the natural cubic spline for the input data points
 * and returns the spline coefficients for each interval between the data points.
 * The spline is defined by the piecewise cubic polynomial:
 *   S(x) = a + b*(x - x_i) + c*(x - x_i)^2 + d*(x - x_i)^3
 * where 'a', 'b', 'c', and 'd' are the coefficients for each segment, and 'x_i'
 * is the x-coordinate of the starting point of the segment.
 */
std::vector<SplineSet> spline(std::vector<double> &x, std::vector<double> &y) {
  // todo: warn if x.size() is less than ???
  size_t n = x.size() - 1;
  std::vector<double> a;
  a.insert(a.begin(), y.begin(), y.end());
  std::vector<double> b(n);
  std::vector<double> d(n);
  std::vector<double> h(n);

  for (size_t i = 0; i < n; ++i) h[i] = (x[i + 1] - x[i]);

  std::vector<double> alpha(n);
  for (size_t i = 1; i < n; ++i) alpha[i] = (3.0 * (a[i + 1] - a[i]) / h[i] - 3.0 * (a[i] - a[i - 1]) / h[i - 1]);

  std::vector<double> c(n + 1);
  std::vector<double> l(n + 1);
  std::vector<double> mu(n + 1);
  std::vector<double> z(n + 1);
  l[0]  = 1.0;
  mu[0] = 0.0;
  z[0]  = 0.0;

  for (size_t i = 1; i < n; ++i) {
    l[i]  = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
  }

  l[n] = 1.0;
  z[n] = 0.0;
  c[n] = 0.0;

  for (size_t j = n - 1; (long)j >= 0; --j) {
    c[j] = z[j] - mu[j] * c[j + 1];
    b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3.0;
    d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
  }

  std::vector<SplineSet> output_set(n + 1);
  for (size_t i = 0; i < n; ++i) {
    output_set[i].a = a[i];
    output_set[i].b = b[i];
    output_set[i].c = c[i];
    output_set[i].d = d[i];
    output_set[i].x = x[i];
  }

  output_set[n].x = x[n]; // the last point is only used to save the last x
  return output_set;
}

void getSlineCurve(std::vector<SplineSet> &cs, std::vector<double> &xsv, std::vector<double> &ysv, int ndiv = 5) {
  for (size_t i = 0; i < cs.size() - 1; ++i) {
    double dx = (cs[i + 1].x - cs[i].x) / (double)ndiv;
    for (double xs = cs[i].x; xs < cs[i + 1].x; xs += dx) {
      double delta = xs - cs[i].x;
      double ys    = cs[i].a + cs[i].b * delta + cs[i].c * delta * delta + cs[i].d * delta * delta * delta;
      xsv.push_back(xs);
      ysv.push_back(ys);
    }
  }
  // last point
  size_t i2    = cs.size() - 1;
  size_t i1    = i2 - 1;
  double delta = cs[i2].x - cs[i1].x;
  double ys    = cs[i1].a + cs[i1].b * delta + cs[i1].c * delta * delta + cs[i1].d * delta * delta * delta;
  xsv.push_back(cs[i2].x);
  ysv.push_back(ys);
}

void derivSpline(std::vector<SplineSet> &cs, std::vector<SplineSet> &csd) {
  SplineSet ss;
  for (size_t i = 0; i < cs.size(); i++) {
    ss.a = cs[i].b;
    ss.b = 2.0 * cs[i].c;
    ss.c = 3.0 * cs[i].d;
    ss.d = 0.0;
    ss.x = cs[i].x;
    csd.push_back(ss);
  }
  /*
  int last = csd.size() - 1;
  csd[last].a = 0.0;
  csd[last].b = 0.0;
  csd[last].c = 0.0;
  csd[last].d = 0.0;
  */
}

} // end of unnamed namespace

/**
@file cubicSpline.hpp
@see
https://en.wikipedia.org/w/index.php?title=Spline_%28mathematics%29&oldid=288288033#Algorithm_for_computing_natural_cubic_splines
Example of usage:
@code{.cpp}
#include <cubicSpline.hpp>
#include <iostream>

int main()
{
        std::vector<double> x(11);
        std::vector<double> y(11);
        for(int i = 0; i < x.size(); ++i) {
                x[i] = i;
                y[i] = sin(i);
        }

        std::vector<SplineSet> cs = spline(x, y);
        for(int i = 0; i < cs.size(); ++i)
                std::cout << cs[i].d << "\t" << cs[i].c << "\t" << cs[i].b << "\t" << cs[i].a << std::endl;
}
@endcode
*/

#if 0

#include <fstream>

int main() {
  std::vector<double> x(8);
  std::vector<double> y(8);
  std::ofstream in("spline_input.txt");
  for (int i = 0; i < x.size(); ++i) {
    x[i] = i*2;
    y[i] = 1.0 / (1.0 + exp(-(x[i] - 8.0)));
    in << x[i] << ' ' << y[i] << '\n';
  }

  std::vector<SplineSet> cs = spline(x, y);

  std::vector<double> xs, ys;
  std::vector<double> xsd, ysd;
  std::vector<SplineSet> csd;
  derivSpline(cs, csd);
  getSlineCurve(csd, xsd, ysd, 10);
  getSlineCurve(cs, xs, ys, 10);
  std::ofstream out("spline_ouput.txt");
  for (int i = 0; i < xs.size(); ++i)
    out << xs[i] << "\t" << ys[i] << "\t" << xsd[i] << "\t" << ysd[i] << std::endl;

  // plot [-1:11] sin(x/2), "test.txt" u 1:2 w p pt 6, cos(x/2), "test.txt" u 3:4 w p pt 6
}

#endif

#endif /* end of include guard: CUBICSPLINE_HPP */

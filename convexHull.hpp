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

#ifndef CONVEXHULL_HPP
#define CONVEXHULL_HPP

#include <algorithm>
#include <vector>

#include "vec2.hpp"

/**
@file convexHull.hpp
@code{.cpp}
#include <convexHull.hpp>
#include <iostream>

int main (int argc, char const *argv[])
{
        std::vector<vec2r> points;
        points.push_back(vec2r(0,0));
        points.push_back(vec2r(1,1));
        points.push_back(vec2r(1,0));
        points.push_back(vec2r(-1,0));

        std::vector<vec2r> ch = convexHull(points);
        for (size_t i = 0 ; i < ch.size() ; i++) {
                std::cout << ch[i] << std::endl;
        }

        return 0;
}
@endcode
*/

namespace {

/// @see    https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
/// @brief  Returns a list of points on the convex hull in counter-clockwise order.
/// @note   The last point in the returned list is the same than the first one.
std::vector<vec2r> convexHull(std::vector<vec2r> &P) {
  int n = P.size(), k = 0;
  if (n == 1) { return P; }
  std::vector<vec2r> H(2 * n);

  // Sort points lexicographically
  std::sort(P.begin(), P.end(), std::less<vec2r>());

  // Build lower hull
  for (int i = 0; i < n; ++i) {
    while (k >= 2 && cross(H[k - 1] - H[k - 2], P[i] - H[k - 2]) <= 0) { k--; }
    H[k++] = P[i];
  }

  // Build upper hull
  for (int i = n - 2, t = k + 1; i >= 0; i--) {
    while (k >= t && cross(H[k - 1] - H[k - 2], P[i] - H[k - 2]) <= 0) { k--; }
    H[k++] = P[i];
  }

  H.resize(k - 1);
  return H;
}

} // end of unnamed namespace

#endif /* end of include guard: CONVEXHULL_HPP */

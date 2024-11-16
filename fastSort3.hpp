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

#ifndef FASTSORT3_HPP
#define FASTSORT3_HPP

#include <utility>

/**
 @brief Inplace-sort the 3 values so that a0 <= a1 <= a2
 @param[in,out] a0 first value
 @param[in,out] a1 second value
 @param[in,out] a2 third value
*/
template <typename T> void fastSort3(T &a0, T &a1, T &a2) {
  if (a0 < a1) {
    if (a1 < a2) {
      return;
    } else {
      if (a0 < a2) {
        std::swap(a1, a2);
      } else {
        std::swap(a0, a2);
        std::swap(a1, a2);
      }
    }
  } else {
    if (a0 < a2) {
      std::swap(a0, a1);
    } else {
      if (a1 < a2) {
        T temp = a0;
        a0 = a1;
        a1 = a2;
        a2 = temp;
      } else {
        std::swap(a0, a2);
      }
    }
  }
}

/**
 @brief Return the indices of the sorted values a0, a1, a2
 @param[in] a0 first value
 @param[in] a1 second value
 @param[in] a2 third value
 @param[out] i0 index of the smallest value
 @param[out] i1 index of the middle value
 @param[out] i2 index of the largest value
*/
template <typename T> void getOrder3(const T a0, const T a1, const T a2, int &i0, int &i1, int &i2) {
  i0 = 0;
  i1 = 1;
  i2 = 2;
  if (a0 < a1) {
    if (a1 < a2) {
      return;
    } else {
      if (a0 < a2) {
        std::swap(i1, i2);
      } else {
        std::swap(i0, i2);
        std::swap(i1, i2);
      }
    }
  } else {
    if (a0 < a2) {
      std::swap(i0, i1);
    } else {
      if (a1 < a2) {
        int temp = i0;
        i0 = i1;
        i1 = i2;
        i2 = temp;
      } else {
        std::swap(i0, i2);
      }
    }
  }
}

#endif /* end of include guard: FASTSORT3_HPP */

#if 0
#include <iostream>
int main(int argc, char const *argv[]) {
  int a = 10, b = 9, c = 9;
  std::cout << "BEFORE SORTING\na = " << a << ", b = " << b << ", c = " << c << '\n';
  int i0, i1, i2;
  getOrder3(a, b, c, i0, i1, i2);
  std::cout << "INDICES\nindex(a) = " << i0 << ", index(b) = " << i1 << ", index(c) = " << i2 << '\n';

  fastSort3<int>(a, b, c);
  std::cout << "AFTER SORTING\na = " << a << ", b = " << b << ", c = " << c << '\n';

  return 0;
}
#endif
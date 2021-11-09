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
 @file
 Usage example:
 @code{.cpp}
 int a = 3, b = 1, c = 2;
 fastSort3<int>(a, b, c);
 @endcode
 The result is: a = 1, b = 2, c = 3
*/

// inplace-sort the 3 values so that a0 <= a1 <= a2
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

#endif /* end of include guard: FASTSORT3_HPP */

#if 0
#include <iostream>
int main(int argc, char const *argv[]) {
  int a = 10, b = 1, c = 2;
  std::cout << "BEFORE SORTING\na = " << a << ", b = " << b << ", c = " << c << '\n';
  fastSort3<int>(a, b, c);
  std::cout << "AFTER SORTING\na = " << a << ", b = " << b << ", c = " << c << '\n';
  return 0;
}
#endif
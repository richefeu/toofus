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

#ifndef COMMON_HPP
#define COMMON_HPP

#include <limits>
#include <string>

// possible other name
// WordSeparator

#define CommBox Common::Instance

typedef double real;
typedef unsigned int uint;

class Common {
public:
  char sep;
  int precision;

  static Common &Instance() {
    static Common inst;
    return inst;
  }

  void set_precision(int p) {
    if (p <= 0)
      precision = std::numeric_limits<double>::digits10 + 1;
    else
      precision = p;
  }

  void sepFromKeyword(std::string &kw) {
    if (kw == "tab")
      sep = '\t';
    else if (kw == "semicolon")
      sep = ';';
    else if (kw == "space")
      sep = ' ';
    else
      sep = ' ';
  }

  std::string keywordFromSep() {
    if (sep == ' ')
      return std::string("space");
    if (sep == '\t')
      return std::string("tab");
    if (sep == ';')
      return std::string("semicolon");
    return std::string("space");
  }

private:
  Common() : sep(' '), precision(std::numeric_limits<double>::digits10 + 1) {}
};

#endif /* end of include guard: COMMON_HPP */

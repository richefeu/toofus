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

#ifndef GRID2_HPP
#define GRID2_HPP

// =============================================
// A two-dimentional grid (Simplified 2D matrix)
// =============================================

// Grid nrows x ncols
template <class T> class grid2 {
public:
  vector<T> Data;
  size_t nrows, ncols;

  grid2() : nrows(0), ncols(0) {}
  grid2(size_t r, size_t c = 1) : Data(r * c, 0), nrows(r), ncols(c) {}

  void resize(size_t r, size_t c = 1) {
    Data.resize(r * c);
    nrows = r;
    ncols = c;
  }

  void fill(const T &val = 0) { std::fill(Data.begin(), Data.end(), val); }

  T &operator()(size_t r, size_t c) { return Data[c * nrows + r]; }

  const T &operator()(size_t r, size_t c) const { return Data[c * nrows + r]; }

  // input/output
  friend ostream &operator<<(ostream &pStr, const grid2 &Grd) {
    for (size_t i = 0; i < Grd.Data.size() - 1; ++i)
      pStr << Grd.Data[i] << ' ';
    pStr << Grd.Data[Grd.Data.size() - 1];
    return pStr;
  }

  friend istream &operator>>(istream &pStr, grid2 &Grd) {
    T val;
    Grd.Data.clear();
    for (size_t i = 0; i < Grd.Data.size() - 1; ++i) {
      pStr >> val;
      Grd.Data.push_back(val);
    }
    return pStr;
  }
};

typedef grid2<double> grid2r;
typedef grid2<int> grid2i;
typedef grid2<unsigned int> grid2ui;
typedef grid2<bool> grid2b;

#endif /* end of include guard: GRID2_HPP */

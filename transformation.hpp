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

#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP

#include "vec3.hpp"

#include <cmath>
#include <iostream>

template <typename T> class Transformation {
public:
  T xx, xy, xz, tx;
  T yx, yy, yz, ty;
  T zx, zy, zz, tz;
  // The last line is always 0 0 0 1

  Transformation() : xx(1), xy(0), xz(0), tx(0), yx(0), yy(1), yz(0), ty(0), zx(0), zy(0), zz(1), tz(0) {}

  Transformation &operator=(const Transformation &M) {
    xx = M.xx;
    xy = M.xy;
    xz = M.xz;
    tx = M.tx;

    yx = M.yx;
    yy = M.yy;
    yz = M.yz;
    ty = M.ty;

    zx = M.zx;
    zy = M.zy;
    zz = M.zz;
    tz = M.tz;
    return (*this);
  }

  Transformation &translate(T x, T y, T z) {
    if (x == (T)0 && y == (T)0 && z == (T)0) {
      return *this;
    }
    Transformation tr;
    tr.tx = x;
    tr.ty = y;
    tr.tz = z;
    *this = *this * tr;
    return *this;
  }

  // Variante matricielle directe
  // https://fr.wikipedia.org/wiki/Formule_d%27Euler%E2%80%93Rodrigues
  Transformation &rotate(T axe_x, T axe_y, T axe_z, T angle) {
    if (angle == (T)0 || (axe_x == (T)0 && axe_y == (T)0 && axe_z == (T)0)) {
      return *this;
    }
    Transformation tr;
    T ca = cos(angle);
    T sa = sin(angle);
    T invn = 1.0f / sqrt(axe_x * axe_x + axe_y * axe_y + axe_z * axe_z);
    T nx = axe_x * invn;
    T ny = axe_y * invn;
    T nz = axe_z * invn;
    T one_ca = (T)1 - ca;

    tr.xx = nx * nx * one_ca + ca;
    tr.xy = nx * ny * one_ca - nz * sa;
    tr.xz = nx * nz * one_ca + ny * sa;

    tr.yx = nx * ny * one_ca + nz * sa;
    tr.yy = ny * ny * one_ca + ca;
    tr.yz = ny * nz * one_ca - nx * sa;

    tr.zx = nx * nz * one_ca - ny * sa;
    tr.zy = ny * nz * one_ca + nx * sa;
    tr.zz = nz * nz * one_ca + ca;

    *this = *this * tr;
    return *this;
  }

  void apply(vec3<T> &v) {
    T vx = xx * v.x + xy * v.y + xz * v.z + tx;
    T vy = yx * v.x + yy * v.y + yz * v.z + ty;
    T vz = zx * v.x + zy * v.y + zz * v.z + tz;
    v.x = vx;
    v.y = vy;
    v.z = vz;
  }

  void reset() {
    xx = 1;
    xy = 0;
    xz = 0;
    tx = 0;
    yx = 0;
    yy = 1;
    yz = 0;
    ty = 0;
    zx = 0;
    zy = 0;
    zz = 1;
    tz = 0;
  }

  void print() {
    std::cout << xx << '\t' << xy << '\t' << xz << '\t' << tx << std::endl;
    std::cout << yx << '\t' << yy << '\t' << yz << '\t' << ty << std::endl;
    std::cout << zx << '\t' << zy << '\t' << zz << '\t' << tz << std::endl;
    std::cout << "0\t0\t0\t1" << std::endl;
  }

  friend Transformation operator*(const Transformation &a, const Transformation &b) {
    Transformation tr;
    tr.xx = a.xx * b.xx + a.xy * b.yx + a.xz * b.zx;
    tr.xy = a.xx * b.xy + a.xy * b.yy + a.xz * b.zy;
    tr.xz = a.xx * b.xz + a.xy * b.yz + a.xz * b.zz;
    tr.tx = a.xx * b.tx + a.xy * b.ty + a.xz * b.tz + a.tx;
    tr.yx = a.yx * b.xx + a.yy * b.yx + a.yz * b.zx;
    tr.yy = a.yx * b.xy + a.yy * b.yy + a.yz * b.zy;
    tr.yz = a.yx * b.xz + a.yy * b.yz + a.yz * b.zz;
    tr.ty = a.yx * b.tx + a.yy * b.ty + a.yz * b.tz + a.ty;
    tr.zx = a.zx * b.xx + a.zy * b.yx + a.zz * b.zx;
    tr.zy = a.zx * b.xy + a.zy * b.yy + a.zz * b.zy;
    tr.zz = a.zx * b.xz + a.zy * b.yz + a.zz * b.zz;
    tr.tz = a.zx * b.tx + a.zy * b.ty + a.zz * b.tz + a.tz;

    return tr;
  }

  // input/output
  friend std::ostream &operator<<(std::ostream &pStr, const Transformation &M) {
    return (pStr << M.xx << ' ' << M.xy << ' ' << M.xz << ' ' << M.tx << ' ' << M.yx << ' ' << M.yy << ' ' << M.yz
                 << ' ' << M.ty << ' ' << M.zx << ' ' << M.zy << ' ' << M.zz << ' ' << M.tz);
  }

  friend std::istream &operator>>(std::istream &pStr, Transformation &M) {
    return (pStr >> M.xx >> M.xy >> M.xz >> M.tx >> M.yx >> M.yy >> M.yz >> M.ty >> M.zx >> M.zy >> M.zz >> M.tz);
  }
};

#endif /* end of include guard: TRANSFORMATION_HPP */

#if 0
#include <iostream>

int main(int argc, char const *argv[]) {

  Transformation<double> T;
  // T.reset();
  //  la derniere tranformation en premier
  T.translate(1., 0., 0.).rotate(0.0, 0.0, 1.0, M_PI / 2).translate(-1., 0., 0.);
  // T.print();
  std::cout << T << std::endl;

  vec3r v(2, 0, 0);
  T.apply(v);
  std::cout << v << std::endl;

  return 0;
}

#endif

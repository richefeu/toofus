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

#ifndef VEC2_HPP
#define VEC2_HPP

// ============================================
// Template class for vectors with 2 components
// ============================================

#include <algorithm>
#include <cmath>
#include <iostream>

/// @brief Vector with 2 components
template <typename T> class vec2 {
public:
  T x, y;

  vec2() : x(0), y(0) {}
  vec2(T X, T Y) : x(X), y(Y) {}
  vec2(const vec2 &v) : x(v.x), y(v.y) {}

  vec2 &operator=(const vec2 &V) {
    x = V.x;
    y = V.y;
    return (*this);
  }

  static vec2 unit_x() { return vec2(1, 0); }
  static vec2 unit_y() { return vec2(0, 1); }
  static vec2 one() { return vec2(1, 1); }
  static vec2 zero() { return vec2(0, 0); }

  void reset() { x = y = 0; }

  /// @brief Set the components of the vector.
  ///
  /// @param X The value to which the x component is set.
  /// @param Y The value to which the y component is set.
  void set(T X, T Y) {
    x = X;
    y = Y;
  }

  /// @brief Set both components of the vector to the given value.
  ///
  /// @param val The value to which both the x and y components are set.
  void set(T val) { x = y = val; }

  /// @brief Check if the vector is null.
  ///
  /// A vector is considered null if both its components are smaller than
  /// the given tolerance.
  ///
  /// @param tol The tolerance value. The default value is 1e-20.
  bool isnull(const T tol = 1e-20) const { return (fabs(x) < tol && fabs(y) < tol); }

  /// @brief Direct access to the vector data.
  ///
  /// This method is very low level and is not recommended for direct use.
  T *c_vec() { return &x; }

  /// @brief Direct access to the vector data.
  ///
  /// This method is very low level and is not recommended for direct use.
  T &operator[](int i) { return *(&x + i); }
  T &operator[](size_t i) { return *(&x + i); }

  /// @brief Direct access to the vector data.
  ///
  /// This method is very low level and is not recommended for direct use.
  const T &operator[](int i) const { return *(&x + i); }
  const T &operator[](size_t i) const { return *(&x + i); }

  // For local frames, the notation n,t and s is more appropriate than x and y
  const T n() const { return x; }
  const T t() const { return y; }

  // Arithmetic operations
  vec2 &operator+=(const vec2 &a) {
    x += a.x;
    y += a.y;
    return *this;
  }

  vec2 &operator-=(const vec2 &a) {
    x -= a.x;
    y -= a.y;
    return *this;
  }

  vec2 &operator*=(T k) {
    x *= k;
    y *= k;
    return *this;
  }

  vec2 &operator/=(T k) {
    x /= k;
    y /= k;
    return *this;
  }

  friend vec2 operator+(const vec2 &a, const vec2 &b) { return vec2(a.x + b.x, a.y + b.y); }

  friend vec2 operator-(const vec2 &a, const vec2 &b) { return vec2(a.x - b.x, a.y - b.y); }

  friend vec2 operator-(const vec2 &a) { return vec2(-a.x, -a.y); }

  friend vec2 operator*(const vec2 &a, T k) { return vec2(a.x * k, a.y * k); }

  friend vec2 operator*(T k, const vec2 &a) { return vec2(a.x * k, a.y * k); }

  friend vec2 operator/(const vec2 &a, T k) { return vec2(a.x / k, a.y / k); }

  // --- Specific external operations ---

  /// Dot product
  friend T operator*(const vec2 &a, const vec2 &b) { return (a.x * b.x + a.y * b.y); }

  /// Multiply each component one another
  friend vec2<T> component_product(const vec2<T> &a, const vec2<T> &b) { return vec2<T>(a.x * b.x, a.y * b.y); }

  /// Find the smallest components
  friend vec2<T> component_min(const vec2<T> &a, const vec2<T> &b) {
    return vec2<T>((a.x < b.x) ? a.x : b.x, (a.y < b.y) ? a.y : b.y);
  }

  /// Find the biggest components
  friend vec2<T> component_max(const vec2<T> &a, const vec2<T> &b) {
    return vec2<T>((a.x > b.x) ? a.x : b.x, (a.y > b.y) ? a.y : b.y);
  }

  /// Absolut value of the components
  friend vec2<T> component_abs(const vec2<T> &a) { return vec2<T>(fabs(a.x), fabs(a.y)); }

  /// Anti-clockwise angle between 2 vectors (from a to b; result in [-pi pi])
  friend T angleBetweenVectors(const vec2<T> &a, const vec2<T> &b) {
    return atan2(a.x * b.y - b.x * a.y, a.x * b.x + a.y * b.y);
  }

  /// Inclination with respect to x (result in [-pi pi])
  friend T inclinationX(const vec2<T> &v) {
    return atan2(v.y, v.x);
  }

  /// angle of vector V with respect to vector Vref (result in [-pi pi])
  /*
  friend T angleBetweenVectors(const vec2<T> & V, const vec2<T> & Vref) {
      T dotp = V * Vref;
      T crossp = cross(V, Vref);
      return atan2(crossp, dotp);
  }
  */

  /// Cross product
  friend T cross(const vec2<T> &a, const vec2<T> &b) { return (a.x * b.y - a.y * b.x); }

  /// Linear interpolation
  friend vec2<T> lerp(double t, const vec2<T> &a, const vec2<T> &b) { return (1.0f - t) * a + t * b; }

  /// Squared length of the vector
  friend T norm2(const vec2 &a) { return a * a; }

  /// Length of the vector
  friend T norm(const vec2 &a) { return sqrt(a * a); }

  T normSup() const { return std::max(std::abs(x), std::abs(y)); }

  T length() const { return norm(*this); }

  /// Normalize and return length (before being normalized)
  T normalize() {
    T N = norm(*this);
    if (N > 0.0)
      *this *= (1.0f / N);
    return N;
  }

  /// Return a normalized vector (without changing 'this' vector)
  vec2 normalized() const {
    vec2 V = *this;
    V.normalize();
    return V;
  }

  /// perform a 90-degree clockwise rotation inplace
  void quarterRightTurn() {
    std::swap(x, y);
    y = -y;
  }
  /// perform a 90-degree anti-clockwise rotation inplace
  void quarterLeftTurn() {
    std::swap(x, y);
    x = -x;
  }

  /// get a 90-degree clockwise rotated vector
  vec2 quarterRightTurned() const { return vec2(y, -x); }

  /// get a 90-degree anti-clockwise rotated vector
  vec2 quarterLeftTurned() const { return vec2(-y, x); }

  /// Determinant
  friend T determinant(const vec2<T> a, const vec2<T> b) { return (a.x * b.y - b.x * a.y); }

  // Comparisons
  bool operator==(const vec2<T> &other) const { return (this->x == other.x && this->y == other.y); }

  bool operator!=(const vec2<T> &other) const { return !(*this == other); }

  // input/output
  friend std::ostream &operator<<(std::ostream &pStr, const vec2 &pV) { return (pStr << pV.x << ' ' << pV.y); }

  friend std::istream &operator>>(std::istream &pStr, vec2 &pV) { return (pStr >> pV.x >> pV.y); }
};

typedef vec2<double> vec2r;
typedef vec2<int> vec2i;
typedef vec2<unsigned int> vec2ui;
typedef vec2<bool> vec2b;

namespace std {
template <class T> struct less<vec2<T>> {
  bool operator()(const vec2<T> &lhs, const vec2<T> &rhs) const {
    if (lhs.x < rhs.x)
      return true;
    else if (lhs.x == rhs.x && lhs.y < rhs.y)
      return true;
    return false;
  }
};
} // namespace std

#endif /* end of include guard: VEC2_HPP */

#if 0
#include <iostream>

int main(int argc, char const *argv[]) {
  vec2r v(1., 2.);
  std::cout << v << '\n';

  //vec2r v2(2., -1.);
  
  //vec2r v2 = v.quarterLeftTurned();
  
  vec2r v2=v;
  v2.quarterRightTurn();
  
  std::cout << v2 << '\n';
  std::cout << cw_angle(v, v2) * 180./M_PI << '\n';

  return 0;
}

#endif

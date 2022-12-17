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

#ifndef MTH_HPP
#define MTH_HPP

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <vector>

namespace Mth {

// All const are with external linkage by default
const double pi = 3.14159265358979323846;
const double invPi = 1.0 / pi;
const double piSqr = pi * pi;
const double pi_2 = pi / 2.0;
const double pi_4 = pi / 4.0;
const double _2pi = 2.0 * pi;
const double _1_3 = 1.0 / 3.0;
const double _4_3 = 4.0 / 3.0;
const double e = 2.71828182845904523536;
const double deg2rad = pi / 180.0;
const double rad2deg = 180.0 / pi;
const double randFactor = 1.0 / (double)(RAND_MAX);

/// @brief Gives angle between 0 and 4,
/// while atan2 gives an angle between -PI and PI
template <typename T> T DiamondAngle(T x, T y) {
  if (y >= 0.0)
    return (x >= 0.0 ? y / (x + y) : 1.0 - x / (-x + y));
  else
    return (x < 0.0 ? 2.0 - y / (-x - y) : 3.0 + x / (x - y));
}

/// @brief Return the sign of a value (1 is positive, -1 is negative)
template <typename T> T sign(T value) { return std::copysign(1, value); }

template <typename T> T map(T value, T minA, T maxA, T minB, T maxB) {
  return (value - minA) / (maxA - minA) * (maxB - minB) + minB;
}

template <typename T> T lerp(T min, T max, T amount) { return min + amount * (max - min); }

template <typename T> T norm(T num, T min, T max) { return (num - min) / (max - min); }

template <typename T> T constrain(T num, T min, T max) {
  if (num < min)
    return min;
  else if (num > max)
    return max;
  return num;
}

template <typename T> T floor0(T value) {
  if (value < 0.0)
    return std::ceil(value);
  else
    return std::floor(value);
}

template <typename T> T round(T value) { return std::floor(value + 0.5); }

template <typename T> T random(T value = 1) { return (std::rand() * randFactor * value); }
template <typename T> T random(T min, T max) { return (min + std::rand() * randFactor * (max - min)); }

template <typename T> T random(std::vector<T> &data) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, data.size() - 1);
  size_t index = static_cast<size_t>(dis(gen));
  return data[index];
}

/// @brief Generate a sobol sequence
/// @param[in] n number of data generated
/// @param[out] x a vector of data of length n
// Note 'inline' is for shut-down the 'unused-function' warning
// and 'static' is to avoid multiple definitions at linkage
inline static void sobolSequence(const size_t n, std::vector<double> &x) {
  const size_t MAXBIT = 30, MAXDIM = 6;
  size_t j, k, l;
  size_t i, im, ipp;
  static size_t mdeg[MAXDIM] = {1, 2, 3, 3, 4, 4};
  static size_t in;
  static std::vector<size_t> ix(MAXDIM);
  static std::vector<size_t *> iu(MAXBIT);
  static size_t ip[MAXDIM] = {0, 1, 1, 2, 1, 4};
  static size_t iv[MAXDIM * MAXBIT] = {1, 1, 1, 1, 1, 1, 3,  1,  3, 3,  1,  1,
                                             5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9};
  static double fac;

  if (n < 0) {
    for (k = 0; k < MAXDIM; k++)
      ix[k] = 0;
    in = 0;
    if (iv[0] != 1)
      return;
    fac = 1.0 / (1 << MAXBIT);
    for (j = 0, k = 0; j < MAXBIT; j++, k += MAXDIM)
      iu[j] = &iv[k];
    for (k = 0; k < MAXDIM; k++) {
      for (j = 0; j < mdeg[k]; j++)
        iu[j][k] <<= (MAXBIT - 1 - j);
      for (j = mdeg[k]; j < MAXBIT; j++) {
        ipp = ip[k];
        i = iu[j - mdeg[k]][k];
        i ^= (i >> mdeg[k]);
        for (l = mdeg[k] - 1; l >= 1; l--) {
          if (ipp & 1)
            i ^= iu[j - l][k];
          ipp >>= 1;
        }
        iu[j][k] = i;
      }
    }
  } else {
    im = in++;
    for (j = 0; j < MAXBIT; j++) {
      if (!(im & 1))
        break;
      im >>= 1;
    }
    if (j >= MAXBIT)
      return; // std::cerr << "MAXBIT too small in sobseq" << std::endl;
    im = j * MAXDIM;
    size_t kmax = (n < MAXDIM) ? n : MAXDIM;
    for (k = 0; k < kmax; k++) {
      ix[k] ^= iv[im + k];
      x[k] = ix[k] * fac;
    }
  }
}

template <typename T> T dist(T x1, T x2) { return std::fabs(x2 - x1); }
template <typename T> T dist(T x1, T y1, T x2, T y2) {
  T dx = x2 - x1;
  T dy = y2 - y1;
  return std::sqrt(dx * dx + dy * dy);
}
template <typename T> T dist(T x1, T y1, T z1, T x2, T y2, T z2) {
  T dx = x2 - x1;
  T dy = y2 - y1;
  T dz = z2 - z1;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}


// QUAKE 3 fast inversed square root
// Thanks to compiler optimisations, using 1.0/sqrt(v) will have
// better performance than these functions
template <typename T> T Q_rsqrt(T number) {
  const T x2 = number * 0.5F;
  const T threehalfs = 1.5F;

  union {
    T f;
    uint32_t i;
  } conv = {.f = number};
  conv.i = 0x5f3759df - (conv.i >> 1);
  conv.f *= threehalfs - (x2 * conv.f * conv.f); // iteration 1
  return conv.f;
}
template <typename T> T Q_accurate_rsqrt(T number) {
  const T x2 = number * 0.5F;
  const T threehalfs = 1.5F;

  union {
    T f;
    uint32_t i;
  } conv = {.f = number};
  conv.i = 0x5f3759df - (conv.i >> 1);
  conv.f *= threehalfs - (x2 * conv.f * conv.f); // iteration 1
  conv.f *= threehalfs - (x2 * conv.f * conv.f); // iteration 2
  return conv.f;
}

/// @brief Compute the mean value and the variance of some data
template <typename T> void MeanAndVariance(std::vector<T> &data, double &mean, double &var) {
  mean = 0.0;
  var = 0.0;
  if (data.size() < 2)
    return;
  for (size_t i = 0; i < data.size(); i++)
    mean += data[i];
  mean /= data.size();
  double s = 0.0;
  for (size_t i = 0; i < data.size(); i++) {
    s = data[i] - mean;
    var += s * s;
  }
  var /= (data.size() - 1);
}

/// @brief Coefficient of variation (CV) or relative standard deviation (RSD)
template <typename T> T RSD(std::vector<T> &data) {
  double mean = 0.0;
  double var = 0.0;
  MeanAndVariance(data, mean, var);
  double CV = 0.0;
  if (fabs(mean) > 1e-20)
    CV = sqrt(var) / mean;
  return CV;
}

} // End of namespace Mth

#if 0
#include <iostream>
int main() {

  /*
  std::cout << "dist = " << Mth::dist(-1.0, -2.0, 6.3 , 2.45) << "\n";
  
  std::vector<double> vec = {2.3,9.3,10.3,98.4,0.34,34.65};
  std::cout << "picked value = " << Mth::random(vec) << "\n";
  std::cout << "picked value = " << Mth::random(vec) << "\n";
  std::cout << "picked value = " << Mth::random(vec) << "\n";
  std::cout << "picked value = " << Mth::random(vec) << "\n";
  */
  
  
  double value = 1234.56876;
  std::cout << 1.0/sqrt(value) << "  vs  " << Mth::Q_rsqrt(value) << "  vs  " << Mth::Q_accurate_rsqrt(value) << '\n';
  
  double a=0.0;
  for (unsigned i = 0 ; i < 1000000000; i++) {
    double v = Mth::random(0.0001,10.);
    double f = 1.0f/sqrt(v);
    //double f = Mth::Q_rsqrt(v);
    a+=f;
  }
  std::cout << "a="<<a<<'\n';



  return 0;
}

#endif

#endif /* end of include guard: MTH_HPP */

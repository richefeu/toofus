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

#ifndef COLORTABLE_HPP
#define COLORTABLE_HPP

#define ALL_BLACK 0
#define VIS5D 1
#define MATLAB_JET 2
#define SAMCET 3
#define RAINBOW 4
#define EMC2000 5
#define BLUE_RED_YELLOW_WHITE 6
#define MATLAB_HOT 7
#define MATLAB_PINK 8
#define GRAYSCALE 9
#define ALL_WHITE 10
#define MATLAB_HSV 11
#define SPECTRUM 12
#define MATLAB_BONE 13
#define MATLAB_SPRING 14
#define MATLAB_SUMMER 15
#define MATLAB_AUTUMN 16
#define MATLAB_WINTER 17
#define MATLAB_COOL 18
#define MATLAB_COPPER 19
#define RANDOM 20

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

struct colorRGBA {
  int r, g, b, a;       // 0 to 255
  float rr, gg, bb, aa; // 0.0 to 1.0 (same data that has been pre-computed)
  colorRGBA() : r(0), g(0), b(0), a(0), rr(0.0), gg(0.0), bb(0.0), aa(0.0) {}

  /**
   * @brief Set the color values.
   *
   * @param R Red component (0 to 255)
   * @param G Green component (0 to 255)
   * @param B Blue component (0 to 255)
   * @param A Alpha component (0 to 255), default is 255
   */
  void set(int R, int G, int B, int A = 255) {
    static const float inv255 = 1.0f / 255.0f;
    r = R;
    g = G;
    b = B;
    a = A;
    rr = (float)r * inv255;
    gg = (float)g * inv255;
    bb = (float)b * inv255;
    aa = (float)a * inv255;
  }

  /**
   * @brief Set the color values using normalized float components.
   *
   * @param RR Red component (0.0 to 1.0)
   * @param GG Green component (0.0 to 1.0)
   * @param BB Blue component (0.0 to 1.0)
   * @param AA Alpha component (0.0 to 1.0), default is 1.0
   *
   * This function updates both integer and float representations of the color
   * values, converting the normalized float components to integer values in
   * the range 0 to 255.
   */
  void set(float RR, float GG, float BB, float AA = 1.0f) {
    rr = RR;
    gg = GG;
    bb = BB;
    aa = AA;
    r = (int)floor(rr * 255);
    g = (int)floor(gg * 255);
    b = (int)floor(bb * 255);
    a = (int)floor(aa * 255);
  }
  
/**
   * @brief Get the color as a hexadecimal string in the format "#RRGGBB".
   *
   * @return A string representing the color in hexadecimal format.
   */
  std::string toHexString() const {
    std::ostringstream oss;
    oss << "#"
        << std::hex << std::setw(2) << std::setfill('0') << r
        << std::hex << std::setw(2) << std::setfill('0') << g
        << std::hex << std::setw(2) << std::setfill('0') << b;
    return oss.str();
  }
};

struct color4ub {
  uint8_t r, g, b, a;
  color4ub() : r(0), g(0), b(0), a(255) {}
};

struct color4f {
  float r, g, b, a;
  color4f() : r(0.0f), g(0.0f), b(0.0f), a(1.0f) {}
};

class ColorTable {
private:
  bool big_endian;
  int size;
  float curvature, bias;
  int rotation;
  bool swap, invert;
  std::vector<unsigned int> table;
  float min, max;
  int tableID;
  float alpha, beta, alphapow;

  /**
   * Pack a color with components in the range 0 to 255 into a single
   * unsigned int, using either big-endian or little-endian byte order.
   *
   * @param R Red component (0 to 255)
   * @param G Green component (0 to 255)
   * @param B Blue component (0 to 255)
   * @param A Alpha component (0 to 255), default is 255
   *
   * @return The packed color value
   */
  unsigned int PACK_COLOR(int R, int G, int B, int A) {
    if (big_endian)
      return ((unsigned int)((R) << 24 | (G) << 16 | (B) << 8 | (A)));
    else
      return ((unsigned int)((A) << 24 | (B) << 16 | (G) << 8 | (R)));
  }

  /**
   * Unpack the red component of a packed color value.
   *
   * @param X Packed color value
   *
   * @return Red component of the color (0 to 255)
   */
  int UNPACK_RED(unsigned int X) {
    if (big_endian)
      return (((X) >> 24) & 0xff);
    else
      return ((X) & 0xff);
  }

  /**
   * Unpack the green component of a packed color value.
   *
   * @param X Packed color value
   *
   * @return Green component of the color (0 to 255)
   */
  int UNPACK_GREEN(unsigned int X) {
    if (big_endian)
      return (((X) >> 16) & 0xff);
    else
      return (((X) >> 8) & 0xff);
  }

  /**
   * Unpack the blue component of a packed color value.
   *
   * @param X Packed color value
   *
   * @return Blue component of the color (0 to 255)
   */
  int UNPACK_BLUE(unsigned int X) {
    if (big_endian)
      return (((X) >> 8) & 0xff);
    else
      return (((X) >> 16) & 0xff);
  }

  /**
   * Unpack the alpha component of a packed color value.
   *
   * @param X Packed color value
   *
   * @return Alpha component of the color (0 to 255)
   */
  int UNPACK_ALPHA(unsigned int X) {
    if (big_endian)
      return ((X) & 0xff);
    else
      return (((X) >> 24) & 0xff);
  }

  /**
   * Compute a gray value from a normalized input value @a s.
   *
   * The mapping is linear from 0.0 to 1.0, and saturated at 0.0 and 1.0.
   *
   * @param s normalized input value
   *
   * @return gray value between 0.0 and 1.0
   */
  double gray(double s) { return s < 0.0 ? 0.0 : (s < 1.0 ? s : 1.0); }

  /**
   * Computes the red component of a 'hot' colormap for a given normalized input value.
   *
   * The function maps an input value 's' to a red intensity value between 0.0 and 1.0.
   * For values of 's' less than 0.0, the output is 0.0. For values of 's' between 0.0 and
   * 3.0/8.0, the output linearly scales from 0.0 to 1.0. For values of 's' greater than or
   * equal to 3.0/8.0, the output is 1.0.
   *
   * @param s Normalized input value (typically between 0.0 and 1.0).
   * @return Red intensity value between 0.0 and 1.0.
   */
  double hot_r(double s) { return s < 0.0 ? 0.0 : (s < 3.0 / 8.0 ? 8.0 / 3.0 * s : 1.0); }

  /**
   * Computes the green component of a 'hot' colormap for a given normalized input value.
   *
   * The function maps an input value 's' to a green intensity value between 0.0 and 1.0.
   * For values of 's' less than 3.0/8.0, the output is 0.0. For values of 's' between
   * 3.0/8.0 and 6.0/8.0, the output linearly scales from 0.0 to 1.0. For values of 's'
   * greater than or equal to 6.0/8.0, the output is 1.0.
   *
   * @param s Normalized input value (typically between 0.0 and 1.0).
   * @return Green intensity value between 0.0 and 1.0.
   */
  double hot_g(double s) { return s < 3.0 / 8.0 ? 0.0 : (s < 6.0 / 8.0 ? 8.0 / 3.0 * (s - 3.0 / 8.0) : 1.0); }

  /**
   * Computes the blue component of a 'hot' colormap for a given normalized input value.
   *
   * The function maps an input value 's' to a blue intensity value between 0.0 and 1.0.
   * For values of 's' less than 6.0/8.0, the output is 0.0. For values of 's' between
   * 6.0/8.0 and 1.0, the output linearly scales from 0.0 to 1.0. For values of 's'
   * greater than or equal to 1.0, the output is 1.0.
   *
   * @param s Normalized input value (typically between 0.0 and 1.0).
   * @return Blue intensity value between 0.0 and 1.0.
   */
  double hot_b(double s) { return s < 6.0 / 8.0 ? 0.0 : (s < 1.0 ? 8.0 / 2.0 * (s - 6.0 / 8.0) : 1.0); }

  /**
   * Compute a cubic polynomial for a given value of x.
   *
   * Given a, b, c, and d, compute a + b*x + c*x^2 + d*x^3.
   *
   * @param a Coefficient of the constant term
   * @param b Coefficient of the linear term
   * @param c Coefficient of the quadratic term
   * @param d Coefficient of the cubic term
   * @param x Value of x
   *
   * @return The computed value of the cubic polynomial
   */
  double cubic(double a, double b, double c, double d, double x) { return a + b * x + c * x * x + d * x * x * x; }

  /**
   * Converts a color specified in HSV (Hue, Saturation, Value) to one in RGB (Red, Green, Blue).
   *
   * @param H Hue value (typically between 0.0 and 6.0).
   * @param S Saturation value (typically between 0.0 and 1.0).
   * @param V Value (brightness) value (typically between 0.0 and 1.0).
   * @param R Output red intensity value between 0.0 and 1.0.
   * @param G Output green intensity value between 0.0 and 1.0.
   * @param B Output blue intensity value between 0.0 and 1.0.
   */
  void HSV_to_RGB(double H, double S, double V, double *R, double *G, double *B) {
    if (S < 5.0e-6) {
      *R = *G = *B = V;
    } else {
      int i = (int)H;
      double f = H - (float)i;
      double p1 = V * (1.0 - S);
      double p2 = V * (1.0 - S * f);
      double p3 = V * (1.0 - S * (1.0 - f));
      switch (i) {
      case 0:
        *R = V;
        *G = p3;
        *B = p1;
        break;
      case 1:
        *R = p2;
        *G = V;
        *B = p1;
        break;
      case 2:
        *R = p1;
        *G = V;
        *B = p3;
        break;
      case 3:
        *R = p1;
        *G = p2;
        *B = V;
        break;
      case 4:
        *R = p3;
        *G = p1;
        *B = V;
        break;
      case 5:
        *R = V;
        *G = p1;
        *B = p2;
        break;
      }
    }
  }

  /**
   * \brief Convert RGB color to HSV
   *
   * This function convert RGB color to HSV color.
   *
   * \param R red component of color
   * \param G green component of color
   * \param B blue component of color
   * \param H pointer to hue value
   * \param S pointer to saturation value
   * \param V pointer to value value
   *
   * \note H value is in range [0, 6.0), S and V values are in range [0, 1.0]
   */
  void RGB_to_HSV(double R, double G, double B, double *H, double *S, double *V) {
    double maxv = R > G ? R : G;
    if (B > maxv)
      maxv = B;
    *V = maxv;
    if (maxv > 0) {
      double minv = R < G ? R : G;
      if (B < minv)
        minv = B;
      *S = 1.0 - double(minv) / maxv;
      if (maxv > minv) {
        if (maxv == R) {
          *H = (G - B) / double(maxv - minv);
          if (*H < 0)
            *H += 6.0;
        } else if (maxv == G)
          *H = 2.0 + (B - R) / double(maxv - minv);
        else
          *H = 4.0 + (R - G) / double(maxv - minv);
      }
    }
  }

public:
  /**
   * \brief Constructor
   *
   * Construct a ColorTable object and initialize internal tables.
   *
   * \param tableID_ table ID
   *
   * Table ID is one of the following:
   * - MATLAB_JET
   * - MATLAB_HSV
   * - MATLAB_HOT
   * - MATLAB_COOL
   * - MATLAB_BONE
   * - MATLAB_WINTER
   * - MATLAB_AUTUMN
   * - MATLAB_SPRING
   * - MATLAB_SUMMER
   * - BINARY
   * - GRAY
   */
  ColorTable(int tableID_ = MATLAB_JET) {
    short int word = 0x0001;
    char *byte = (char *)&word;
    big_endian = (byte[0] ? false : true);

    rotation = 0.0;
    bias = 0.0;
    curvature = 0.0;
    min = 0.0;
    max = 1.0;
    size = 256;
    swap = false;
    tableID = tableID_;
    invert = false;
    alphapow = 0.0;
    alpha = 0.0;
    beta = 0.0;

    Rebuild();
  }

  /**
   * \brief Rebuilds the color table based on the current configuration.
   *
   * This function clears and reconstructs the internal color table using
   * the specified size, rotation, bias, curvature, and other parameters.
   * The colors are generated based on the selected color table ID, which
   * determines the color mapping scheme (e.g., MATLAB_JET, MATLAB_HSV, etc.).
   * The function handles various color schemes including black, white, grayscale,
   * and several MATLAB named schemes, among others. It also applies adjustments
   * based on inversion, alpha blending, and gamma correction using the beta parameter.
   */
  void Rebuild() {
    double s, t, gamma;
    int r, g, b, a;

    if (!table.empty())
      table.clear();
    table.reserve(size);

    for (int i = 0; i < size; i++) {

      if (size > 1) {
        if (i + rotation < 0)
          s = (double)(i + rotation + size) / (double)(size - 1);
        else if (i + rotation > size - 1)
          s = (double)(i + rotation - size) / (double)(size - 1);
        else
          s = (double)(i + rotation) / (double)(size - 1);
      } else
        s = 0.0;

      if (swap)
        s = 1.0 - s;

      switch (tableID) {
      case 0: // all black
        r = g = b = 0;
        break;
      case 1: // vis5d
        t = (curvature + 1.4) * (s - (1.0 + bias) / 2.0);
        r = (int)(128.0 + 127.0 * atan(7.0 * t) / 1.57);
        g = (int)(128.0 + 127.0 * (2 * exp(-7 * t * t) - 1));
        b = (int)(128.0 + 127.0 * atan(-7.0 * t) / 1.57);
        break;
      case 2: { // matlab "jet"
        double ii = (double)(s - bias) * 128.0;
        if (ii < 0)
          ii = 0;
        if (ii > 128)
          ii = 128;
        double rr = ii <= 46 ? 0.0 : ii >= 111 ? -0.03125 * (ii - 111) + 1.0 : ii >= 78 ? 1. : 0.03125 * (ii - 46);
        double gg = ii <= 14 || ii >= 111 ? 0.
                    : ii >= 79            ? -0.03125 * (ii - 111)
                    : ii <= 46            ? 0.03125 * (ii - 14)
                                          : 1.0;
        double bb = ii >= 79 ? 0.0 : ii >= 47 ? -0.03125 * (ii - 79) : ii <= 14 ? 0.03125 * (ii - 14) + 1.0 : 1.0;
        r = (int)(rr * 255.0);
        g = (int)(gg * 255.0);
        b = (int)(bb * 255.0);
      } break;
      case 3: // lucie, samcef (?)
        if (s - bias <= 0.0) {
          r = 0;
          g = 0;
          b = 255;
        } else if (s - bias <= 0.40) {
          r = 0;
          g = (int)((s - bias) * 637.5);
          b = (int)(255. - (s - bias) * 637.5);
        } else if (s - bias <= 0.60) {
          r = (int)(1275.0 * (s - bias - 0.4));
          g = 255;
          b = 0;
        } else if (s - bias <= 1.0) {
          r = 255;
          g = (int)(255.0 - 637.5 * (s - bias - 0.6));
          b = 0;
        } else {
          r = 255;
          g = 0;
          b = 0;
        }
        break;
      case 4: // rainbow
        if (s - bias <= 0.0) {
          r = 0;
          g = 0;
          b = 255;
        } else if (s - bias <= 0.25 + curvature) {
          curvature = (curvature == -0.25f) ? -0.26f : curvature;
          r = 0;
          g = (int)((s - bias) * (255.0 / (0.25 + curvature)));
          b = 255;
        } else if (s - bias <= 0.50) {
          curvature = (curvature == 0.25f) ? 0.26f : curvature;
          r = 0;
          g = 255;
          b = (int)(255. - (255. / (0.25 - curvature)) * (s - bias - 0.25 - curvature));
        } else if (s - bias <= 0.75 - curvature) {
          curvature = (curvature == 0.25f) ? 0.26f : curvature;
          r = (int)((s - bias - 0.5) * (255.0 / (0.25 - curvature)));
          g = 255;
          b = 0;
        } else if (s - bias <= 1.) {
          curvature = (curvature == -0.25f) ? -0.26f : curvature;
          r = 255;
          g = (int)(255.0 - (255.0 / (0.25 + curvature)) * (s - bias - 0.75 + curvature));
          b = 0;
        } else {
          r = 255;
          g = 0;
          b = 0;
        }
        break;
      case 5: // emc2000 (rainbow with black and white)
        if (s - bias <= 0.) {
          r = 0;
          g = 0;
          b = 0;
        } else if (s - bias <= 0.2) {
          r = (int)(57 * (1 - 100 * ((s - bias) - 0.1) * ((s - bias) - 0.1)));
          g = 0;
          b = (int)((s - bias) * (255. / 0.2));
        } else if (s - bias <= 0.3624) {
          r = 0;
          g = (int)((s - bias - 0.2) * (255. / 0.1624));
          b = 255;
        } else if (s - bias <= 0.50) {
          r = 0;
          g = 255;
          b = (int)(255. - (255. / 0.1376) * (s - bias - 0.3624));
        } else if (s - bias <= 0.6376) {
          r = (int)((s - bias - 0.5) * (255. / 0.1376));
          g = 255;
          b = 0;
        } else if (s - bias <= 0.8) {
          r = 255;
          g = (int)(255. - (255. / 0.1624) * (s - bias - 0.6376));
          b = 0;
        } else if (s - bias <= 1.0) {
          r = 255;
          g = (int)((255. / 0.2) * (s - bias - 0.8));
          b = (int)(-3187.66 * (s - bias) * (s - bias) + 7012.76 * (s - bias) - 3570.61);
        } else {
          r = 255;
          g = 255;
          b = 255;
        }
        break;
      case 6: // darkblue->red->yellow->white
        r = (int)(255.0 * cubic(-0.0506169, 2.81633, -1.87033, 0.0524573, s - bias));
        g = (int)(255.0 * cubic(0.0485868, -1.26109, 6.3074, -4.12498, s - bias));
        b = (int)(255.0 * cubic(0.364662, 1.50814, -7.36756, 6.51847, s - bias));
        break;
      case 7: // matlab "hot"
        r = (int)(255.0 * hot_r(s - bias));
        g = (int)(255.0 * hot_g(s - bias));
        b = (int)(255.0 * hot_b(s - bias));
        break;
      case 8: // matlab "pink"
        r = (int)(255.0 * sqrt((2.0 * gray(s - bias) + hot_r(s - bias)) / 3.0));
        g = (int)(255.0 * sqrt((2.0 * gray(s - bias) + hot_g(s - bias)) / 3.0));
        b = (int)(255.0 * sqrt((2.0 * gray(s - bias) + hot_b(s - bias)) / 3.0));
        break;
      case 9: // grayscale
        if (s - bias <= 0.0) {
          r = g = b = 0;
        } else if (s - bias <= 1.) {
          r = g = b = (int)(255 * (1.0 - curvature) * (s - bias));
        } else {
          r = g = b = (int)(255 * (1.0 - curvature));
        }
        break;
      case 10: // all white
        r = g = b = 255;
        break;
      case 11: { // matlab "hsv"
        double H = 6.0 * s + 1.0e-10, R = 0.0, G = 0.0, B = 0.0;
        HSV_to_RGB(H, 1.0, 1.0, &R, &G, &B);
        r = (int)(255 * R);
        g = (int)(255 * G);
        b = (int)(255 * B);
      } break;
      case 12: { // spectrum (truncated hsv)
        double H = 5.0 * s + 1.e-10, R = 0.0, G = 0.0, B = 0.0;
        HSV_to_RGB(H, 1.0, 1.0, &R, &G, &B);
        r = (int)(255 * R);
        g = (int)(255 * G);
        b = (int)(255 * B);
      } break;
      case 13: // matlab "bone"
        r = (int)(255.0 * (7.0 * gray(s - bias) + hot_b(s - bias)) / 8.0);
        g = (int)(255.0 * (7.0 * gray(s - bias) + hot_g(s - bias)) / 8.0);
        b = (int)(255.0 * (7.0 * gray(s - bias) + hot_r(s - bias)) / 8.0);
        break;
      case 14: // matlab "spring"
        r = (int)(255.0 * 1.0);
        g = (int)(255.0 * gray(s - bias));
        b = (int)(255.0 * (1.0 - gray(s - bias)));
        break;
      case 15: // matlab "summer"
        r = (int)(255.0 * gray(s - bias));
        g = (int)(255.0 * (0.5 + gray(s - bias) / 2.0));
        b = (int)(255.0 * 0.4);
        break;
      case 16: // matlab "autumn"
        r = (int)(255.0 * 1.0);
        g = (int)(255.0 * gray(s - bias));
        b = (int)(255.0 * 0.0);
        break;
      case 17: // matlab "winter"
        r = (int)(255.0 * 0.0);
        g = (int)(255.0 * gray(s - bias));
        b = (int)(255.0 * (0.5 + (1.0 - gray(s - bias)) / 2.0));
        break;
      case 18: // matlab "cool"
        r = (int)(255.0 * gray(s - bias));
        g = (int)(255.0 * (1.0 - gray(s - bias)));
        b = (int)(255.0 * 1.0);
        break;
      case 19: // matlab "copper"
        r = (int)(255.0 * MIN(1.0, gray(s - bias) * 1.25));
        g = (int)(255.0 * MIN(1.0, gray(s - bias) * 0.7812));
        b = (int)(255.0 * MIN(1.0, gray(s - bias) * 0.4975));
        break;
      case 20: // Random colors
        r = (int)(rand() / (double)RAND_MAX * 255.0);
        g = (int)(rand() / (double)RAND_MAX * 255.0);
        b = (int)(rand() / (double)RAND_MAX * 255.0);
        break;
      default:
        r = g = b = 0;
        break;
      }

      float aa = 1.0;
      if (alphapow)
        aa = static_cast<float>(pow((float)(s ? s : 1.0e-10), (float)alphapow));
      a = (int)(255.0 * aa * alpha);

      if (beta) {
        if (beta > 0.0)
          gamma = 1.0 - beta;
        else
          gamma = 1.0 / (1.001 + beta); // beta is thresholded to [-1,1]
        r = (int)(255.0 * pow((float)r / 255.0, gamma));
        g = (int)(255.0 * pow((float)g / 255.0, gamma));
        b = (int)(255.0 * pow((float)b / 255.0, gamma));
      }

      if (invert) {
        r = 255 - r;
        g = 255 - g;
        b = 255 - b;
      }

      // clamp to [0,255]
      r = r < 0 ? 0 : (r > 255 ? 255 : r);
      g = g < 0 ? 0 : (g > 255 ? 255 : g);
      b = b < 0 ? 0 : (b > 255 ? 255 : b);
      a = a < 0 ? 0 : (a > 255 ? 255 : a);

      table[i] = PACK_COLOR(r, g, b, a);
    }
  }

  /**
   * Rebuilds the color table by interpolating RGBA values.
   *
   * This function clears the current color table and rebuilds it by interpolating
   * between the specified colors at given positions. The interpolation is done
   * linearly based on the positions in `cpos` and the colors in `cols`.
   *
   * @param cpos A vector of integer positions representing key points in the color table.
   *             These positions must be within the valid range [0, size) and must be sorted
   *             in strictly increasing order.
   * @param cols A vector of vectors, where each inner vector contains four integers
   *             representing the RGBA values of a color. The size of `cols` must match
   *             the size of `cpos`, and each RGBA value should be in the range [0, 255].
   *
   * The function checks for various conditions, such as:
   * - `cpos` and `cols` must have the same size.
   * - Elements in `cpos` must be within the valid range and sorted.
   *
   * If any condition is violated, an error message is printed and the function returns
   * without modifying the table.
   */
  void rebuild_interp_rgba(std::vector<int> cpos, std::vector<std::vector<int>> cols) {
    if (!table.empty())
      table.clear();
    table.reserve(size);

    if (cpos.size() != cols.size()) {
      std::cerr << "In colorTable::rebuild, cpos and cols are not of the same size\n";
      return;
    }

    for (size_t i = 0; i < cpos.size(); i++) {
      if (cpos[i] < 0 || cpos[i] >= size) {
        std::cerr << "In colorTable::rebuild, cpos element out of range\n";
        return;
      }
    }

    for (size_t i = 1; i < cpos.size(); i++) {
      if (cpos[i] <= cpos[i - 1]) {
        std::cerr << "In colorTable::rebuild, cpos not sorted\n";
        return;
      }

      colorRGBA col;
      for (int c = cpos[i - 1]; c <= cpos[i]; c++) {
        float w = (float)(c - cpos[i - 1]) / (float)(cpos[i] - cpos[i - 1]);
        col.r = (int)((1.0 - w) * (float)cols[i - 1][0] + w * (float)cols[i][0]);
        col.g = (int)((1.0 - w) * (float)cols[i - 1][1] + w * (float)cols[i][1]);
        col.b = (int)((1.0 - w) * (float)cols[i - 1][2] + w * (float)cols[i][2]);
        col.a = (int)((1.0 - w) * (float)cols[i - 1][3] + w * (float)cols[i][3]);

        // clamp to [0,255]
        col.r = col.r < 0 ? 0 : (col.r > 255 ? 255 : col.r);
        col.g = col.g < 0 ? 0 : (col.g > 255 ? 255 : col.g);
        col.b = col.b < 0 ? 0 : (col.b > 255 ? 255 : col.b);
        col.a = col.a < 0 ? 0 : (col.a > 255 ? 255 : col.a);

        table[c] = PACK_COLOR(col.r, col.g, col.b, col.a);
      }
    }
  }

  /**
   * \brief Rebuild color table with linear interpolation in HSV color space
   *
   * \param cpos Vector of color positions
   * \param cols Vector of color values at positions
   *
   * The color table is cleared and rebuilt with linear interpolation in HSV
   * color space. If the size of \a cpos and \a cols are not the same, or if a
   * position in \a cpos is out of range, the function prints an error message
   * and returns without modifying the color table. If the positions in \a cpos
   * are not sorted, the function prints an error message and returns without
   * modifying the color table.
   */
  void rebuild_interp_hsv(std::vector<int> cpos, std::vector<std::vector<int>> cols) {
    if (!table.empty())
      table.clear();
    table.reserve(size);

    if (cpos.size() != cols.size()) {
      std::cerr << "In colorTable::rebuild, cpos and cols are not of the same size\n";
      return;
    }

    for (size_t i = 0; i < cpos.size(); i++) {
      if (cpos[i] < 0 || cpos[i] >= size) {
        std::cerr << "In colorTable::rebuild, cpos element out of range\n";
        return;
      }
    }

    for (size_t i = 1; i < cpos.size(); i++) {
      if (cpos[i] <= cpos[i - 1]) {
        std::cerr << "In colorTable::rebuild, cpos not sorted\n";
        return;
      }

      colorRGBA col;
      for (int c = cpos[i - 1]; c <= cpos[i]; c++) {
        float w = (float)(c - cpos[i - 1]) / (float)(cpos[i] - cpos[i - 1]);

        double R1 = cols[i - 1][0] / 255.0;
        double G1 = cols[i - 1][1] / 255.0;
        double B1 = cols[i - 1][2] / 255.0;
        double R2 = cols[i][0] / 255.0;
        double G2 = cols[i][1] / 255.0;
        double B2 = cols[i][2] / 255.0;
        double H1 = 0.0, S1 = 0.0, V1 = 0.0, H2 = 0.0, S2 = 0.0, V2 = 0.0;
        RGB_to_HSV(R1, G1, B1, &H1, &S1, &V1);
        RGB_to_HSV(R2, G2, B2, &H2, &S2, &V2);
        double H = (1.0 - w) * H1 + w * H2;
        double S = (1.0 - w) * S1 + w * S2;
        double V = (1.0 - w) * V1 + w * V2;
        double R = 0.0, G = 0.0, B = 0.0;
        HSV_to_RGB(H, S, V, &R, &G, &B);
        col.r = (int)(255 * R);
        col.g = (int)(255 * G);
        col.b = (int)(255 * B);

        // clamp to [0,255]
        col.r = col.r < 0 ? 0 : (col.r > 255 ? 255 : col.r);
        col.g = col.g < 0 ? 0 : (col.g > 255 ? 255 : col.g);
        col.b = col.b < 0 ? 0 : (col.b > 255 ? 255 : col.b);
        col.a = 255;

        table[c] = PACK_COLOR(col.r, col.g, col.b, col.a);
      }
    }
  }

  /**
   * Save the color table to a PPM file. The color table is displayed vertically,
   * with the lowest color index at the top and the highest at the bottom. The
   * width of the image is determined as the smallest power of 2 that is larger
   * than the height divided by 15.
   *
   * @param name the name of the PPM file
   */
  void savePpm(const char *name) {
    int h = size;
    int w = (int)(h / 15 + 1);

    std::ofstream file(name, std::ios::binary);
    file << "P6" << '\n';
    file << w << " " << h << '\n';
    file << "255" << '\n';

    unsigned char r, g, b;
    int icol;
    for (int y = 0; y < h; y++) {
      icol = h - y - 1;
      for (int x = 0; x < w; x++) {
        r = (unsigned char)UNPACK_RED(table[icol]);
        g = (unsigned char)UNPACK_GREEN(table[icol]);
        b = (unsigned char)UNPACK_BLUE(table[icol]);
        file.write((const char *)&r, sizeof(unsigned char));
        file.write((const char *)&g, sizeof(unsigned char));
        file.write((const char *)&b, sizeof(unsigned char));
      }
    }
  }

  void setTableID(int id) { tableID = id; }
  void setSwap(bool s) { swap = s; }
  void setInvert(bool i) { invert = i; }

  /**
   * Set the minimum and maximum values for the color table. These values are
   * used to determine the color for each value in the table. The color table
   * is evenly divided between the minimum and maximum values.
   *
   * @param Min the minimum value
   * @param Max the maximum value
   */
  void setMinMax(float Min, float Max) {
    min = Min;
    max = Max;
  }
  void setSize(int Size) { size = Size; }
  void setBias(float Bias) { bias = Bias; }
  void setCurvature(float Curv) { curvature = Curv; }
  void setRotation(int Rot) { rotation = Rot; }

  int getSize() { return size; }
  float getMin() { return min; }
  float getMax() { return max; }

  /**
   * Retrieves the RGBA color corresponding to the given value from the color table.
   *
   * This function maps a float value to its respective color in the color table
   * by normalizing the value based on the minimum and maximum range. The resulting
   * color is stored in the provided colorRGBA structure.
   *
   * @param value The float value to map to a color, expected to be within the
   *              range [min, max].
   * @param col A pointer to a colorRGBA structure where the resulting color will
   *            be stored. The structure's fields `r`, `g`, `b`, and `a` are set
   *            to the corresponding red, green, blue, and alpha components of the color.
   *
   * The function ensures that the color index is clamped between the bounds of
   * the color table and extracts the color components using predefined unpacking
   * functions.
   */
  void getRGB(float value, colorRGBA *col) {
    unsigned int i;

    float pos = (value - min) / (max - min);                // in ]0.0 1.0[
    i = (unsigned int)(floor(pos * (float)(size - 2))) + 1; // in [1 size-2]

    if (value <= min)
      i = 0;
    else if (value >= max)
      i = size - 1;

    col->r = UNPACK_RED(table[i]);
    col->g = UNPACK_GREEN(table[i]);
    col->b = UNPACK_BLUE(table[i]);
    col->a = UNPACK_ALPHA(table[i]);
  }

  /**
   * Retrieves the 8-bit RGBA color corresponding to the given value from the color table.
   *
   * This function maps a float value to its respective color in the color table
   * by normalizing the value based on the minimum and maximum range. The resulting
   * color is stored in the provided color4ub structure as 8-bit unsigned integers.
   *
   * @param value The float value to map to a color, expected to be within the
   *              range [min, max].
   * @param col A pointer to a color4ub structure where the resulting color will
   *            be stored. The structure's fields `r`, `g`, `b`, and `a` are set
   *            to the corresponding red, green, blue, and alpha components of the color.
   *
   * The function ensures that the color index is clamped between the bounds of
   * the color table and extracts the color components using predefined unpacking
   * functions.
   */
  void getColor4ub(float value, color4ub *col) {
    unsigned int i;

    float pos = (value - min) / (max - min);                // in ]0.0 1.0[
    i = (unsigned int)(floor(pos * (float)(size - 2))) + 1; // in [1 size-2]

    if (value <= min)
      i = 0;
    else if (value >= max)
      i = size - 1;

    col->r = static_cast<uint8_t>(UNPACK_RED(table[i]));
    col->g = static_cast<uint8_t>(UNPACK_GREEN(table[i]));
    col->b = static_cast<uint8_t>(UNPACK_BLUE(table[i]));
    col->a = static_cast<uint8_t>(UNPACK_ALPHA(table[i]));
  }

  /**
   * Retrieves the floating-point RGBA color corresponding to the given value from the color table.
   *
   * This function maps a float value to its respective color in the color table
   * by normalizing the value based on the minimum and maximum range. The resulting
   * color is stored in the provided color4f structure as floats in the range [0.0, 1.0].
   *
   * @param value The float value to map to a color, expected to be within the
   *              range [min, max].
   * @param col A pointer to a color4f structure where the resulting color will
   *            be stored. The structure's fields `r`, `g`, `b`, and `a` are set
   *            to the corresponding red, green, blue, and alpha components of the color.
   *
   * The function ensures that the color index is clamped between the bounds of
   * the color table and extracts the color components using predefined unpacking
   * functions.
   */
  void getColor4f(float value, color4f *col) {
    static const float inv255 = 1.0f / 255.0f;
    unsigned int i;

    float pos = (value - min) / (max - min);                // in ]0.0 1.0[
    i = (unsigned int)(floor(pos * (float)(size - 2))) + 1; // in [1 size-2]

    if (value <= min)
      i = 0;
    else if (value >= max)
      i = size - 1;

    col->r = static_cast<float>(inv255 * (float)UNPACK_RED(table[i]));
    col->g = static_cast<float>(inv255 * (float)UNPACK_GREEN(table[i]));
    col->b = static_cast<float>(inv255 * (float)UNPACK_BLUE(table[i]));
    col->a = static_cast<float>(inv255 * (float)UNPACK_ALPHA(table[i]));
  }

  /**
   * Retrieves a random color from a set of 8 predefined colors.
   *
   * The colors are evenly spaced in the RGB cube, and are well
   * distinguishable from each other.
   *
   * @param col A pointer to a colorRGBA structure where the resulting color will
   *            be stored. The structure's fields `r`, `g`, `b`, and `a` are set
   *            to the corresponding red, green, blue, and alpha components of the color.
   */
  void getRandomRGB8(colorRGBA *col) {
    static const std::vector<unsigned int> tb = {
        PACK_COLOR(225, 205, 90, 255), PACK_COLOR(170, 114, 160, 255), PACK_COLOR(230, 116, 98, 255),
        PACK_COLOR(94, 129, 148, 255), PACK_COLOR(226, 221, 147, 255), PACK_COLOR(221, 175, 185, 255),
        PACK_COLOR(208, 143, 72, 255), PACK_COLOR(118, 169, 155, 255),
    };
    unsigned int i = (unsigned int)floor(rand() / (double)RAND_MAX * 8);
    if (i >= 8)
      i = 7;
    col->r = UNPACK_RED(tb[i]);
    col->g = UNPACK_GREEN(tb[i]);
    col->b = UNPACK_BLUE(tb[i]);
    col->a = UNPACK_ALPHA(tb[i]);
  }

  /**
   * Retrieves a cyclic color from a set of 8 predefined colors.
   *
   * The colors are evenly spaced in the RGB cube, and are well
   * distinguishable from each other.
   *
   * @param col A pointer to a colorRGBA structure where the resulting color will
   *            be stored. The structure's fields `r`, `g`, `b`, and `a` are set
   *            to the corresponding red, green, blue, and alpha components of the color.
   *
   * The color index is incremented on each call and wraps around to 0
   * when it reaches 8.
   */
  void getCyclicRGB8(colorRGBA *col) {
    static const std::vector<unsigned int> tb = {
        PACK_COLOR(225, 205, 90, 255), PACK_COLOR(170, 114, 160, 255), PACK_COLOR(230, 116, 98, 255),
        PACK_COLOR(94, 129, 148, 255), PACK_COLOR(226, 221, 147, 255), PACK_COLOR(221, 175, 185, 255),
        PACK_COLOR(208, 143, 72, 255), PACK_COLOR(118, 169, 155, 255),
    };
    static unsigned int count = 0;
    col->r = UNPACK_RED(tb[count]);
    col->g = UNPACK_GREEN(tb[count]);
    col->b = UNPACK_BLUE(tb[count]);
    col->a = UNPACK_ALPHA(tb[count]);
    count++;
    if (count >= 8)
      count = 0;
  }

  /**
   * Prints the contents of the color table to the standard output.
   *
   * For each color in the table, the color index and its red, green,
   * blue, and alpha components are printed to the standard output.
   *
   * The purpose of this function is to provide a simple way to verify
   * the correctness of the color table. It is not intended to be used
   * in performance-critical code.
   */
  void Print() {
    int i, r, g, b, a;

    for (i = 0; i < size; i++) {
      r = UNPACK_RED(table[i]);
      g = UNPACK_GREEN(table[i]);
      b = UNPACK_BLUE(table[i]);
      a = UNPACK_ALPHA(table[i]);
      std::cout << i << " [" << r << ", " << g << ", " << b << ", " << a << "]" << std::endl;
    }
  }
};

#endif /* end of include guard: COLORTABLE_HPP */

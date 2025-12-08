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

#ifndef MAT4_HPP
#define MAT4_HPP

#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>

#include "vec2.hpp"

template <typename T> class mat4 {

public:
  T xx, xy;
  T yx, yy;

  mat4() : xx(0), xy(0), yx(0), yy(0) {}
  // mat4(mat4sym & m): xx(m.xx), xy(m.xy), yx(m.xy), yy(m.yy) { }
  mat4(const T XX, const T XY, const T YX, const T YY) : xx(XX), xy(XY), yx(YX), yy(YY) {}
  mat4(const T M[]) : xx(M[0]), xy(M[1]), yx(M[2]), yy(M[3]) {}
  mat4(const mat4 &M) : xx(M.xx), xy(M.xy), yx(M.yx), yy(M.yy) {}

  mat4 &operator=(const mat4 &M) {
    xx = M.xx;
    xy = M.xy;
    yx = M.yx;
    yy = M.yy;
    return (*this);
  }

  // Constants
  static mat4 unit() {
    return mat4(1, 0, 0, 1);
  }
  static mat4 zero() {
    return mat4(1, 1, 1, 1);
  }
  static mat4 one() {
    return mat4(1, 1, 1, 1);
  }

  /// Sets all elements of the matrix to 0.
  void reset() {
    xx = xy = yx = yy = 0;
  }

  void reset(const double val) {
    xx = xy = val;
    yx = yy = val;
  }

  /// Sets the diagonal elements of the matrix.
  ///
  /// \param[in] XX First element of the diagonal.
  /// \param[in] YY Second element of the diagonal.
  void set_diag(const double XX, const double YY) {
    xx = XX;
    yy = YY;
  }

  T &operator[](int i) {
    return *(&xx + i);
  }

  T &operator[](size_t i) {
    return *(&xx + i);
  }

  const T &operator[](int i) const {
    return *(&xx + i);
  }
  const T &operator[](size_t i) const {
    return *(&xx + i);
  }

  T &at(int line, int column) {
    return *(&xx + 2 * line + column);
  }
  const T &at(int line, int column) const {
    return *(&xx + 2 * line + column);
  }

  T *c_mtx() {
    return &xx;
  }

  /// Sets all elements of the matrix to zero.
  void setZero() {
    xx = xy = yx = yy = 0.0;
  }

  /// Sets the matrix to the identity matrix.
  ///
  /// This method modifies the matrix such that the diagonal elements
  /// are set to 1 and the off-diagonal elements are set to 0, effectively
  /// transforming it into an identity matrix.
  void setIdentity() {
    xx = yy = 1.0;
    xy = yx = 0.0;
  }

  /// Returns a new matrix which is the transpose of the current matrix.
  ///
  /// \return A new matrix which is the transpose of the current matrix.
  mat4 transposed() {
    return mat4(xx, yx, xy, yy);
  }

  /// Sets the current matrix to its transpose.
  ///
  /// This method modifies the matrix by swapping the off-diagonal elements,
  /// effectively transforming it into its transpose.
  void transpose() {
    std::swap(xy, yx);
  }

  /// only for symmetric matrix
  void eigenvalues(double &v1, double &v2, bool &swapped) const {
    v1 = 0.5 * (xx + yy) + sqrt((0.5 * (xx - yy)) * (0.5 * (xx - yy)) + xy * xy);
    v2 = 0.5 * (xx + yy) - sqrt((0.5 * (xx - yy)) * (0.5 * (xx - yy)) + xy * xy);
    if (v2 > v1) {
      double swap = v1;
      v1          = v2;
      v2          = swap;
      swapped     = true;
    }
  }

  /// Computes the eigenvalues and eigenvectors of the matrix and stores them in
  /// the given matrices \p V and \p D. The eigenvalues are stored in \p D in
  /// ascending order, and the corresponding eigenvectors are stored as columns
  /// in \p V. The eigenvectors are normalized.
  ///
  /// \param[out] V Matrix to store the eigenvectors.
  /// \param[out] D Matrix to store the eigenvalues.
  void eigen(mat4 &V, mat4 &D) {
    double TT  = xx + yy;
    double det = xx * yy - xy * yx;
    double L1  = 0.5 * TT + sqrt(0.25 * TT * TT - det); // eigenval
    double L2  = 0.5 * TT - sqrt(0.25 * TT * TT - det); // eigenval
    D.xx       = L1;
    D.xy = D.yx = 0.0;
    D.yy        = L2;
    // Eigenvectors organized vertically
    if (yx != 0) {
      V.xx = L1 - yy;
      V.yx = yx;
      V.xy = L2 - yy;
      V.yy = yx;
    } else if (xy != 0) {
      V.xx = xy;
      V.yx = L1 - xx;
      V.xy = xy;
      V.yy = L2 - xx;
    } else if (xy == 0 and yx == 0) {
      V.xx = 1;
      V.yx = 0;
      V.xy = 0;
      V.yy = 1;
    }

    // Normalizing vectors
    vec2r v1(V.xx, V.yx);
    vec2r v2(V.xy, V.yy);
    v1 = v1 / sqrt(v1.x * v1.x + v1.y * v1.y); // use normalized()...
    v2 = v2 / sqrt(v2.x * v2.x + v2.y * v2.y);

    // Putting them back in the V matrix
    V.reset();
    V.xx = v1.x;
    V.yx = v1.y;
    V.xy = v2.x;
    V.yy = v2.y;
  }

  /**
   * @brief Compute eigenvectors (stored as columns in V) and corresponding eigenvalues (D) by assuming the matrix is
   * double and symmetric
   *
   * See section 11.1 of Numerical Recipes in C for more information.
   *
   * @param V matrix to store eigenvectors as columns
   * @param D matrix to store eigenvalues
   * @return the number of rotations, or -1 if the matrix is not symmetric
   */
  int sym_eigen(mat4 &V, mat4 &D) const {
    int rot = 0;
    vec2r B;
    vec2r Z;

    // Save the input matrix in orig, use new matrix inp
    mat4 A = *this;
    // Set vectors to the identity matrix
    V.xx = 1;
    V.xy = 0;
    V.yx = 0;
    V.yy = 1;
    // Set B and D values to the diagonal of the input matrix
    B.x = D.xx = A.xx;
    B.y = D.yy = A.yy;

    // Rotate until off-diagonal elements of input matrix are zero
    for (int sweep = 0; sweep++ < 50;) {
      double sum = fabs(A.xy);
      double thresh;

      if (fabs(sum) < 1.0e-15) return rot;

      thresh   = (sweep < 4) ? sum * 0.2 / 4.0 : 0.0; // First three sweeps?
      double g = 100.0 * fabs(A.xy);                  // TBC!!

      // After 4 sweeps, skip the rotation if the
      // off-diagonal element is small.
      if ((sweep > 4) && (g < 1.0e-15)) A.xy = 0.0;
      else if (fabs(A.xy) > thresh) {

        double h = D.yy - D.xx;
        double c, s, t; // cosine, sine, tangent of rotation angle
        double tau;

        if (g < 1.0e-20) t = A.xy / h;
        else {
          double theta = 0.5 * h / A.xy;
          t            = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
          if (theta < 0.0) t = -t;
        }

        c   = 1.0 / sqrt(1.0 + t * t); // cosine of rotation angle
        s   = t * c;                   // sine of rotation angle
        tau = s / (1.0 + c);

        h = t * A.xy;
        Z.x -= h;
        Z.y += h;
        D.xx -= h;
        D.yy += h;
        A.xy = 0.0;

        g    = V.xx;
        h    = V.xy;
        V.xx = g - s * (h + g * tau);
        V.xy = h + s * (g - h * tau);

        g    = V.yx;
        h    = V.yy;
        V.yx = g - s * (h + g * tau);
        V.yy = h + s * (g - h * tau);

        rot++;
      }

      // Set the eigen values
      B += Z;
      D.xx = B.x;
      D.yy = B.y;
      // D = B;
      Z.x = 0.0;
      Z.y = 0.0;
    }
    return -1; // Non-normal return - too many rotations
  }

  bool inverse() {
    double det = xx * yy - xy * yx;
    if (fabs(det) < 1.0e-20) return false; // inverse cannot be calculated

    double swap = xx;
    xx          = yy;
    yy          = swap;

    double inv_det = 1.0 / det;
    xx *= inv_det;
    xy *= -inv_det;
    yy *= inv_det;
    yx *= -inv_det;

    return true;
  }

  mat4 get_inverse() {
    double det = xx * yy - xy * yx;
    // if (fabs(det) < 1.0e-20) return false; // inverse cannot be calculated
    double xx1(xx), xy1(xy), yx1(yx), yy1(yy);
    double swap = xx1;
    xx1         = yy1;
    yy1         = swap;

    double inv_det = 1.0 / det;
    xx1 *= inv_det;
    xy1 *= -inv_det;
    yy1 *= inv_det;
    yx1 *= -inv_det;

    return mat4(xx1, xy1, yx1, yy1);
  }

  T normSup() const {
    return std::max({std::abs(xx), std::abs(xy), std::abs(yx), std::abs(yy)});
  }

  T det() const {
    return (xx * yy - xy * yx);
  }

  T trace() const {
    return (xx + yy);
  }

  void svd(mat4 &U, mat4 &S, mat4 &V) const {
    // taken from http://www.lucidarme.me/?p=4802
    // U matrix
    double val1 = xx * yx + xy * yy;
    double val2 = xx * xx + xy * xy - yx * yx - yy * yy;
    double val3 = xx * xy + yx * yy;
    double val4 = xx * xx - xy * xy + yx * yx - yy * yy;

    double theta = 0.5 * atan2(2 * val1, val2);
    U.xx         = cos(theta);
    U.xy         = -sin(theta);
    U.yx         = sin(theta);
    U.yy         = cos(theta);

    // Singular value matrix (S)
    double S1 = xx * xx + xy * xy + yx * yx + yy * yy;
    double S2 = sqrt(val2 * val2 + 4.0 * val1 * val1);
    // singular values
    double sv1 = sqrt((S1 + S2) / 2.0);
    double sv2 = sqrt((S1 - S2) / 2.0);

    S.xx = sv1;
    S.yy = sv2;
    S.xy = 0;
    S.yx = 0;

    // V matrix
    double phi = 0.5 * atan2(2.0 * val3, val4);
    double s11 = (xx * cos(theta) + yx * sin(theta)) * cos(phi) + (xy * cos(theta) + yy * sin(theta)) * sin(phi);
    double s22 = (xx * sin(theta) - yx * cos(theta)) * sin(phi) + (-xy * sin(theta) + yy * cos(theta)) * cos(phi);
    V.xx       = s11 / fabs(s11) * cos(phi);
    V.xy       = -s22 / fabs(s22) * sin(phi);
    V.yx       = s11 / fabs(s11) * sin(phi);
    V.yy       = s22 / fabs(s22) * cos(phi);
  }

  bool square_root(mat4 &SqR) const {
    double tau   = xx + yy;
    double delta = xx * yy - xy * yx;

    if (delta == 0.0) return false;

    double s = sqrt(delta);
    double t = sqrt(tau + 2 * s);

    SqR.xx = (xx + s) / t;
    SqR.xy = (xy) / t;
    SqR.yx = (yx) / t;
    SqR.yy = (yy + s) / t;

    return true;
  }

  // =======================
  //  Arithmetic operations
  // =======================

  mat4 &operator+=(const mat4 &a) {
    xx += a.xx;
    xy += a.xy;
    yx += a.yx;
    yy += a.yy;
    return *this;
  }

  mat4 &operator-=(const mat4 &a) {
    xx -= a.xx;
    xy -= a.xy;
    yx -= a.yx;
    yy -= a.yy;
    return *this;
  }

  mat4 &operator*=(double k) {
    xx *= k;
    xy *= k;
    yx *= k;
    yy *= k;
    return *this;
  }

  mat4 &operator/=(double k) {
    xx /= k;
    xy /= k;
    yx /= k;
    yy /= k;
    return *this;
  }

  // Comparisons
  bool operator==(const mat4 &other) const {
    return (this->xx == other.xx && this->xy == other.xy && this->yx == other.yx && this->yy == other.yy);
  }

  bool operator!=(const mat4 &other) const {
    return !(*this == other);
  }

  // =========
  //  FRIENDS
  // =========

  friend mat4 operator+(const mat4 &a, const mat4 &b) {
    return mat4(a.xx + b.xx, a.xy + b.xy, a.yx + b.yx, a.yy + b.yy);
  }

  friend mat4 operator-(const mat4 &a, const mat4 &b) {
    return mat4(a.xx - b.xx, a.xy - b.xy, a.yx - b.yx, a.yy - b.yy);
  }

  friend mat4 operator-(const mat4 &a) {
    return mat4(-a.xx, -a.xy, -a.yx, -a.yy);
  }

  friend mat4 operator*(const mat4 &a, double k) {
    return mat4(k * a.xx, k * a.xy, k * a.yx, k * a.yy);
  }

  friend mat4 operator*(double k, const mat4 &a) {
    return mat4(k * a.xx, k * a.xy, k * a.yx, k * a.yy);
  }

  friend mat4 operator/(const mat4 &a, double k) {
    return mat4(a.xx / k, a.xy / k, a.yx / k, a.yy / k);
  }

  friend vec2r operator*(const mat4 &a, const vec2r &b) {
    return vec2r(a.xx * b.x + a.xy * b.y, a.yx * b.x + a.yy * b.y);
  }

  friend vec2r operator*(const vec2r &b, const mat4 &a) {
    return vec2r(a.xx * b.x + a.xy * b.y, a.yx * b.x + a.yy * b.y);
  }

  friend mat4 operator*(const mat4 &a, const mat4 &b) {
    return mat4(a.xx * b.xx + a.xy * b.yx, a.xx * b.xy + a.xy * b.yy, a.yx * b.xx + a.yy * b.yx,
                a.yx * b.xy + a.yy * b.yy);
  }

  friend std::ostream &operator<<(std::ostream &pStr, const mat4 &pV) {
    return (pStr << pV.xx << ' ' << pV.xy << ' ' << pV.yx << ' ' << pV.yy);
  }

  friend std::istream &operator>>(std::istream &pStr, mat4 &M) {
    return (pStr >> M.xx >> M.xy >> M.yx >> M.yy);
  }

  // --- Style flags ---  
  static const int NoOption{0};
  static const int ColoredBrackets{1 << 0};
  static const int WithSeparators{1 << 1};
  static const int Compact{1 << 2};
  static const int Scientific{1 << 3};

  void fancyPrint(int opts = NoOption, int precision = 3) const {
    const std::string reset = "\033[0m";
    const std::string blue  = "\033[34m";

    bool coloredBrackets = (opts & ColoredBrackets);
    bool withSeparators  = (opts & WithSeparators);
    bool compact         = (opts & Compact);
    bool scientific      = (opts & Scientific);

    std::string sep = withSeparators ? (compact ? "|" : " | ") : (compact ? " " : "   ");

    // --- Step 1: format all values as strings ---
    auto fmt = [&](double v) {
      std::ostringstream oss;
      oss << std::setprecision(precision);
      if (scientific) oss << std::scientific;
      else oss << std::fixed;
      oss << v;
      return oss.str();
    };

    std::string s_xx = fmt(xx);
    std::string s_xy = fmt(xy);
    std::string s_yx = fmt(yx);
    std::string s_yy = fmt(yy);

    // --- Step 2: compute column widths ---
    size_t width_col1 = std::max(s_xx.size(), s_yx.size());
    size_t width_col2 = std::max(s_xy.size(), s_yy.size());

    // --- Step 3: select big brackets ---
    std::string topLeft     = "⎡";
    std::string bottomLeft  = "⎣";
    std::string topRight    = "⎤";
    std::string bottomRight = "⎦";

    if (coloredBrackets) {
      topLeft     = blue + topLeft + reset;
      bottomLeft  = blue + bottomLeft + reset;
      topRight    = blue + topRight + reset;
      bottomRight = blue + bottomRight + reset;
    }

    // --- Step 4: print rows aligned with one big bracket ---
    std::cout << std::fixed << std::setprecision(precision);

    std::ostringstream row1, row2;
    row1 << " " << std::setw(static_cast<int>(width_col1)) << s_xx << sep << std::setw(static_cast<int>(width_col2))
         << s_xy << " ";
    row2 << " " << std::setw(static_cast<int>(width_col1)) << s_yx << sep << std::setw(static_cast<int>(width_col2))
         << s_yy << " ";

    std::cout << topLeft << row1.str() << topRight << "\n";
    std::cout << bottomLeft << row2.str() << bottomRight << "\n";
  }
};

// predefined typedefs
typedef mat4<double> mat4r;
typedef mat4<float> mat4f;
typedef mat4<int> mat4i;
typedef mat4<unsigned int> mat4ui;
typedef mat4<bool> mat4b;

#endif /* end of include guard: MAT4_HPP */
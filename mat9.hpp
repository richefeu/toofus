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

#ifndef MAT9_HPP
#define MAT9_HPP

/// @file
/// @brief 3 by 3 matrix
/// @author Vincent Richefeu <Vincent.Richefeu@3sr-grenoble.fr>,
/// Lab 3SR, Grenoble University

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>

#include "vec3.hpp"

/// Matrix 3x3
template <typename T> class mat9 {
public:
  T xx, xy, xz;
  T yx, yy, yz;
  T zx, zy, zz;

  /**
   * @brief Default constructor for a 3x3 matrix.
   *
   * Initializes all elements of the matrix to zero, creating a zero matrix.
   */
  mat9() : xx(0), xy(0), xz(0), yx(0), yy(0), yz(0), zx(0), zy(0), zz(0) {}

  /**
   * @brief Constructs a 3x3 matrix from its elements.
   *
   * This constructor is useful when the elements of the matrix are known
   * at compile-time.
   *
   * @param XX The element in the first row, first column.
   * @param XY The element in the first row, second column.
   * @param XZ The element in the first row, third column.
   * @param YX The element in the second row, first column.
   * @param YY The element in the second row, second column.
   * @param YZ The element in the second row, third column.
   * @param ZX The element in the third row, first column.
   * @param ZY The element in the third row, second column.
   * @param ZZ The element in the third row, third column.
   */
  mat9(const T XX, const T XY, const T XZ, const T YX, const T YY, const T YZ, const T ZX, const T ZY, const T ZZ)
      : xx(XX), xy(XY), xz(XZ), yx(YX), yy(YY), yz(YZ), zx(ZX), zy(ZY), zz(ZZ) {}

  /**
   * @brief Constructs a 3x3 matrix with all elements set to a specified value.
   *
   * This constructor is useful for initializing a matrix with a single
   * value, such as zero or one.
   *
   * @param val The value to which all elements of the matrix are set.
   */
  mat9(const T val) : xx(val), xy(val), xz(val), yx(val), yy(val), yz(val), zx(val), zy(val), zz(val) {}

  /**
   * @brief Constructs a 3x3 matrix from three column vectors.
   *
   * The columns of the matrix are set to the values of the three
   * provided vectors.
   *
   * @param col1 The first column of the matrix.
   * @param col2 The second column of the matrix.
   * @param col3 The third column of the matrix.
   */
  mat9(const vec3<T> &col1, const vec3<T> &col2, const vec3<T> &col3) {
    set(col1, col2, col3);
  }

  /**
   * @brief Constructs a 3x3 matrix from an array.
   *
   * This constructor initializes the matrix elements using the values
   * provided in a single-dimensional array. The array is expected to have
   * at least 9 elements, where each element corresponds to a specific
   * position in the 3x3 matrix, filled in row-major order.
   *
   * @param M An array of at least 9 elements representing the matrix values.
   */
  mat9(const T M[]) {
    xx = M[0];
    xy = M[1];
    xz = M[2];
    yx = M[3];
    yy = M[4];
    yz = M[5];
    zx = M[6];
    zy = M[7];
    zz = M[8];
  }

  /**
   * @brief Copy constructor.
   *
   * This constructor takes a mat9 as an argument and creates a new mat9
   * that is a copy of the input matrix.
   *
   * @param M The matrix to be copied.
   */
  mat9(const mat9 &M) : xx(M.xx), xy(M.xy), xz(M.xz), yx(M.yx), yy(M.yy), yz(M.yz), zx(M.zx), zy(M.zy), zz(M.zz) {}

  mat9 &operator=(const mat9 &M) {
    xx = M.xx;
    xy = M.xy;
    xz = M.xz;
    yx = M.yx;
    yy = M.yy;
    yz = M.yz;
    zx = M.zx;
    zy = M.zy;
    zz = M.zz;
    return (*this);
  }

  /**
   * @brief Returns a 3x3 zero matrix.
   *
   * This static function creates a mat9 object representing a 3x3 matrix
   * with all elements initialized to zero.
   *
   * @return A 3x3 zero matrix.
   */
  static mat9 zero() {
    return mat9(0, 0, 0, 0, 0, 0, 0, 0, 0);
  }

  /**
   * @brief Returns a 3x3 identity matrix.
   *
   * This static function creates a mat9 object representing a 3x3 matrix
   * with all diagonal elements set to one and all other elements set to
   * zero, i.e. the identity matrix.
   *
   * @return The 3x3 identity matrix.
   */
  static mat9 unit() {
    return mat9(1, 0, 0, 0, 1, 0, 0, 0, 1);
  }

  /**
   * @brief Returns a 3x3 matrix with all elements set to one.
   *
   * This static function creates a mat9 object representing a 3x3 matrix
   * with all elements set to one.
   *
   * @return A 3x3 matrix with all elements set to one.
   */
  static mat9 one() {
    return mat9(1, 1, 1, 1, 1, 1, 1, 1, 1);
  }

  /**
   * Sets the elements of the matrix to the specified values.
   *
   * @param XX The element in the first row, first column.
   * @param XY The element in the first row, second column.
   * @param XZ The element in the first row, third column.
   * @param YX The element in the second row, first column.
   * @param YY The element in the second row, second column.
   * @param YZ The element in the second row, third column.
   * @param ZX The element in the third row, first column.
   * @param ZY The element in the third row, second column.
   * @param ZZ The element in the third row, third column.
   */
  void set(const T XX, const T XY, const T XZ, const T YX, const T YY, const T YZ, const T ZX, const T ZY, const T ZZ) {
    xx = XX;
    xy = XY;
    xz = XZ;
    yx = YX;
    yy = YY;
    yz = YZ;
    zx = ZX;
    zy = ZY;
    zz = ZZ;
  }

  /**
   * @brief Sets the matrix columns using three vec3 vectors.
   *
   * This method takes three vectors representing the columns of the matrix
   * and assigns their components to the corresponding elements of the matrix.
   *
   * @param col1 The first column vector.
   * @param col2 The second column vector.
   * @param col3 The third column vector.
   */
  void set(const vec3<T> &col1, const vec3<T> &col2, const vec3<T> &col3) {
    xx = col1.x;
    xy = col2.x;
    xz = col3.x;
    yx = col1.y;
    yy = col2.y;
    yz = col3.y;
    zx = col1.z;
    zy = col2.z;
    zz = col3.z;
  }

  /**
   * @brief Sets the specified column of the matrix with the given vector.
   *
   * This method takes an index and a vector and assigns the vector's components
   * to the specified column of the matrix. The index determines which column
   * is set to the vector's values.
   *
   * @param i The index of the column to be set (0-based).
   * @param col The vector whose components are assigned to the column.
   */
  void set_col(int i, const vec3<T> &col) {
    *(&xx + i)     = col.x;
    *(&xx + 3 + i) = col.y;
    *(&xx + 6 + i) = col.z;
  }

  /**
   * @brief Resets the matrix to zero.
   *
   * This method sets all the elements of the matrix to zero.
   */
  void reset() {
    xx = xy = xz = 0.0;
    yx = yy = yz = 0.0;
    zx = zy = zz = 0.0;
  }

  /**
   * @brief Resets the matrix to a specified value.
   *
   * This method sets all the elements of the matrix to the provided value.
   *
   * @param val The value to which all elements of the matrix are set.
   */
  void reset(const T val) {
    xx = xy = xz = val;
    yx = yy = yz = val;
    zx = zy = zz = val;
  }

  /**
   * @brief Sets the diagonal elements of the matrix.
   *
   * This method takes three parameters and assigns them to the
   * diagonal elements of the matrix. The diagonal elements of the
   * matrix are the elements where the row and column indices are
   * equal.
   *
   * @param XX The element in the first row, first column.
   * @param YY The element in the second row, second column.
   * @param ZZ The element in the third row, third column.
   */
  void set_diag(const T XX, const T YY, const T ZZ) {
    xx = XX;
    yy = YY;
    zz = ZZ;
  }

  /**
   * @brief Accesses the element at the specified index.
   *
   * This is an overloaded operator which allows the user to access
   * elements of the matrix with the square bracket operator. The index
   * argument is the same as the index of the element in a 1D array of
   * size 9.
   *
   * @param i The index of the element to access (0-based).
   * @return The element at the specified index.
   */
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

  /**
   * @brief Accesses the element at the specified row and column.
   *
   * This method takes two parameters and returns the element of the
   * matrix at the specified row and column.
   *
   * @param line The row of the element to access (0-based).
   * @param column The column of the element to access (0-based).
   * @return The element at the specified row and column.
   */
  T &at(int line, int column) {
    return *(&xx + 3 * line + column);
  }
  const T &at(int line, int column) const {
    return *(&xx + 3 * line + column);
  }

  /**
   * @brief Returns a pointer to the first element of the matrix.
   *
   * This method returns a pointer to the first element of the matrix.
   * This is useful when working with functions that require a pointer
   * to the first element of a matrix.
   *
   * @return A pointer to the first element of the matrix.
   */
  T *c_mtx() {
    return &xx;
  }

  /**
   * @brief Adds two 3x3 matrices element-wise.
   *
   * This function takes two matrices as input and returns a new matrix
   * where each element is the sum of the corresponding elements in the
   * input matrices.
   *
   * @param a The first matrix operand.
   * @param b The second matrix operand.
   * @return A new matrix resulting from the element-wise addition of matrices a and b.
   */
  friend mat9 operator+(const mat9 &a, const mat9 &b) {
    return mat9(a.xx + b.xx, a.xy + b.xy, a.xz + b.xz, a.yx + b.yx, a.yy + b.yy, a.yz + b.yz, a.zx + b.zx, a.zy + b.zy,
                a.zz + b.zz);
  }

  /**
   * @brief Subtracts two 3x3 matrices element-wise.
   *
   * This function takes two matrices as input and returns a new matrix
   * where each element is the difference of the corresponding elements
   * in the input matrices.
   *
   * @param a The minuend matrix.
   * @param b The subtrahend matrix.
   * @return A new matrix resulting from the element-wise subtraction of matrices a and b.
   */
  friend mat9 operator-(const mat9 &a, const mat9 &b) {
    return mat9(a.xx - b.xx, a.xy - b.xy, a.xz - b.xz, a.yx - b.yx, a.yy - b.yy, a.yz - b.yz, a.zx - b.zx, a.zy - b.zy,
                a.zz - b.zz);
  }

  /**
   * @brief Returns the negation of a 3x3 matrix.
   *
   * This function takes a matrix as input and returns a new matrix
   * where each element is the negation of the corresponding element in the
   * input matrix.
   *
   * @param a The matrix to be negated.
   * @return A new matrix resulting from the negation of matrix a.
   */
  friend mat9 operator-(const mat9 &a) {
    return mat9(-a.xx, -a.xy, -a.xz, -a.yx, -a.yy, -a.yz, -a.zx, -a.zy, -a.zz);
  }

  /**
   * @brief Scales a 3x3 matrix by a scalar.
   *
   * This function takes a matrix and a scalar as input and returns a new
   * matrix where each element is the product of the corresponding element
   * in the input matrix and the scalar.
   *
   * @param a The matrix to be scaled.
   * @param k The scalar factor.
   * @return A new matrix resulting from the scaling of matrix a.
   */
  friend mat9 operator*(const mat9 &a, T k) {
    return mat9(k * a.xx, k * a.xy, k * a.xz, k * a.yx, k * a.yy, k * a.yz, k * a.zx, k * a.zy, k * a.zz);
  }

  /**
   * @brief Returns the product of two 3x3 matrices.
   *
   * This function takes two matrices as input and returns a new
   * matrix where each element is the dot product of the corresponding
   * row in the first matrix and column in the second matrix.
   *
   * @param a The first matrix.
   * @param b The second matrix.
   * @return A new matrix resulting from the multiplication of matrices a and b.
   */
  friend mat9 operator*(const mat9 &a, const mat9 &b) {
    return mat9(a.xx * b.xx + a.xy * b.yx + a.xz * b.zx, a.xx * b.xy + a.xy * b.yy + a.xz * b.zy,
                a.xx * b.xz + a.xy * b.yz + a.xz * b.zz, a.yx * b.xx + a.yy * b.yx + a.yz * b.zx,
                a.yx * b.xy + a.yy * b.yy + a.yz * b.zy, a.yx * b.xz + a.yy * b.yz + a.yz * b.zz,
                a.zx * b.xx + a.zy * b.yx + a.zz * b.zx, a.zx * b.xy + a.zy * b.yy + a.zz * b.zy,
                a.zx * b.xz + a.zy * b.yz + a.zz * b.zz);
  }

  /**
   * @brief Computes the inner product between two 3x3 matrices.
   *
   * This function takes two matrices as input and returns a scalar value
   * which is the sum of the products of the corresponding elements of the
   * input matrices.
   *
   * The inner product is computed as the trace of the product of the first
   * matrix with the transpose of the second matrix.
   *
   * @param a The first matrix.
   * @param b The second matrix.
   * @return The inner product of matrices a and b.
   */
  friend double inner_product(const mat9 &a, const mat9 &b) {
    mat9 bt = b;
    bt.transpose();
    return (a * bt).trace();
  }

  /**
   * @brief Computes the Hadamard product of two 3x3 matrices.
   *
   * This function takes two matrices as input and returns a new matrix
   * where each element is the product of the corresponding elements of
   * the input matrices.
   *
   * @param a The first matrix.
   * @param b The second matrix.
   * @return A new matrix resulting from the Hadamard product of matrices a and b.
   */
  friend mat9 hadamard_product(const mat9 &a, const mat9 &b) {
    return mat9(a.xx * b.xx, a.xy * b.xy, a.xz * b.xz, a.yx * b.yx, a.yy * b.yy, a.yz * b.yz, a.zx * b.zx, a.zy * b.zy,
                a.zz * b.zz);
  }

  /**
   * @brief Computes the spherical part of a 3x3 matrix.
   *
   * This function takes a matrix as input and returns a new matrix
   * where all diagonal elements are equal to one-third of the trace of
   * the input matrix, and all off-diagonal elements are zero.
   *
   * @param a The input matrix.
   * @return A new matrix representing the spherical part of the input matrix.
   */
  friend mat9 spheric(const mat9 &a) {
    double tr_3 = (a.trace() / 3.0);
    return mat9(tr_3, 0, 0, 0, tr_3, 0, 0, 0, tr_3);
  }

  /**
   * @brief Computes the deviatoric part of a 3x3 matrix.
   *
   * This function takes a matrix as input and returns a new matrix
   * where all diagonal elements are equal to the corresponding diagonal
   * elements of the input matrix minus one-third of the trace of the input
   * matrix, and all off-diagonal elements are equal to the corresponding
   * off-diagonal elements of the input matrix.
   *
   * @param a The input matrix.
   * @return A new matrix representing the deviatoric part of the input matrix.
   */
  friend mat9 deviatoric(const mat9 &a) {
    return a - spheric(a);
  }

  /**
   * @brief Computes the product of a scalar and a 3x3 matrix.
   *
   * This function takes a scalar and a matrix as input and returns a new matrix
   * where each element is the product of the scalar and the corresponding element
   * of the input matrix.
   *
   * @param k The scalar.
   * @param a The input matrix.
   * @return A new matrix resulting from the product of the scalar and the input matrix.
   */
  friend mat9 operator*(T k, const mat9 &a) {
    return mat9(k * a.xx, k * a.xy, k * a.xz, k * a.yx, k * a.yy, k * a.yz, k * a.zx, k * a.zy, k * a.zz);
  }

  /**
   * @brief Computes the product of a 3x3 matrix and a 3-element vector.
   *
   * This function takes a matrix and a vector as input and returns a new vector
   * where each element is the dot product of the corresponding row of the matrix
   * and the input vector.
   *
   * @param a The input matrix.
   * @param v The input vector.
   * @return A new vector resulting from the product of the matrix and the input vector.
   */
  friend vec3<T> operator*(const mat9 &a, const vec3<T> &v) {
    return vec3<T>(a.xx * v.x + a.xy * v.y + a.xz * v.z, a.yx * v.x + a.yy * v.y + a.yz * v.z,
                   a.zx * v.x + a.zy * v.y + a.zz * v.z);
  }

  /**
   * @brief Computes the division of a 3x3 matrix and a scalar.
   *
   * This function takes a matrix and a scalar as input and returns a new matrix
   * where each element is the quotient of the corresponding element of the
   * input matrix and the scalar.
   *
   * @param a The input matrix.
   * @param K The scalar.
   * @return A new matrix resulting from the division of the input matrix and the scalar.
   */
  friend mat9 operator/(const mat9 &a, T K) {
    T k = 0.0;
    if (K != 0.0) k = 1.0 / K;
    return mat9(k * a.xx, k * a.xy, k * a.xz, k * a.yx, k * a.yy, k * a.yz, k * a.zx, k * a.zy, k * a.zz);
  }

  /**
   * @brief Adds the elements of another matrix to this matrix.
   *
   * This operator modifies the current matrix by adding the corresponding
   * elements of the input matrix to it.
   *
   * @param a The matrix to be added to the current matrix.
   */
  void operator+=(const mat9 &a) {
    xx += a.xx;
    xy += a.xy;
    xz += a.xz;
    yx += a.yx;
    yy += a.yy;
    yz += a.yz;
    zx += a.zx;
    zy += a.zy;
    zz += a.zz;
  }

  /**
   * @brief Subtracts the elements of another matrix from this matrix.
   *
   * This operator modifies the current matrix by subtracting the corresponding
   * elements of the input matrix from it.
   *
   * @param a The matrix to be subtracted from the current matrix.
   */
  void operator-=(const mat9 &a) {
    xx -= a.xx;
    xy -= a.xy;
    xz -= a.xz;
    yx -= a.yx;
    yy -= a.yy;
    yz -= a.yz;
    zx -= a.zx;
    zy -= a.zy;
    zz -= a.zz;
  }

  /**
   * @brief Multiplies the elements of this matrix by a scalar.
   *
   * This operator modifies the current matrix by multiplying each of its
   * elements by the scalar input.
   *
   * @param k The scalar to multiply the elements of this matrix with.
   */
  void operator*=(T k) {
    xx *= k;
    xy *= k;
    xz *= k;
    yx *= k;
    yy *= k;
    yz *= k;
    zx *= k;
    zy *= k;
    zz *= k;
  }

  /**
   * @brief Divides the elements of this matrix by a scalar.
   *
   * This operator modifies the current matrix by dividing each of its
   * elements by the scalar input, if the scalar is non-zero. If the scalar
   * is zero, the division is not performed to avoid division by zero.
   *
   * @param K The scalar to divide the elements of this matrix by.
   */
  void operator/=(T K) {
    T k = 0.0;
    if (K != 0.0) { k = 1.0 / K; }
    xx *= k;
    xy *= k;
    xz *= k;
    yx *= k;
    yy *= k;
    yz *= k;
    zx *= k;
    zy *= k;
    zz *= k;
  }

  /**
   * @brief Sets all elements of the matrix to zero.
   *
   * This function modifies the current matrix by setting each of its elements
   * to zero, effectively resetting it to a zero matrix.
   */
  void setZero() {
    xx = xy = xz = yx = yy = yz = zx = zy = zz = 0.0;
  }

  /**
   * @brief Sets the matrix to the identity matrix.
   *
   * This function modifies the current matrix by setting each of its elements
   * to their identity matrix values, effectively resetting it to the identity
   * matrix.
   */
  void setIdentity() {
    xx = yy = zz = 1.0;
    xy = xz = yx = yz = zx = zy = 0.0;
  }

  /**
   * @brief Symmetrizes the matrix.
   *
   * This function modifies the current matrix by adding the corresponding
   * elements of its transpose to it, effectively symmetrizing it.
   */
  void symmetrize() {
    xy = 0.5 * (xy + yx);
    yx = xy;
    xz = 0.5 * (xz + zx);
    zx = xz;
    yz = 0.5 * (yz + zy);
    zy = yz;
  }

  /**
   * @brief Computes the supremum norm of this matrix.
   *
   * This function computes the supremum norm of this matrix, which is the
   * maximum of the absolute values of its elements.
   *
   * @return The supremum norm of this matrix.
   */
  T normSup() const {
    return std::max({std::abs(xx), std::abs(xy), std::abs(xz), std::abs(yx), std::abs(yy), std::abs(yz), std::abs(zx),
                     std::abs(zy), std::abs(zz)});
  }

  /**
   * @brief Computes the determinant of this matrix.
   *
   * This function computes the determinant of this matrix by the standard
   * formula.
   *
   * @return The determinant of this matrix.
   */
  T det() const {
    return (xx * (yy * zz - zy * yz) - yx * (xy * zz - zy * xz) + zx * (xy * yz - yy * xz));
  }

  /**
   * @brief Computes the trace of this matrix.
   *
   * The trace of a matrix is the sum of its diagonal elements.
   *
   * @return The trace of this matrix.
   */
  T trace() const {
    return (xx + yy + zz);
  }

  /**
   * @brief Computes the transpose of this matrix in-place.
   *
   * This function simply swaps the off-diagonal elements of this matrix to
   * compute the transpose in-place.
   */
  void transpose() {
    std::swap(xy, yx);
    std::swap(xz, zx);
    std::swap(yz, zy);
  }

  /**
   * @brief Computes the transpose of this matrix.
   *
   * This function returns a new matrix containing the transpose of this matrix.
   * The transpose of a matrix is the matrix whose rows are the columns of the
   * original matrix and whose columns are the rows of the original matrix.
   *
   * @return The transpose of this matrix.
   */
  mat9<T> transposed() const {
    return mat9<T>(xx, yx, zx, xy, yy, zy, xz, yz, zz);
  }

  /**
   * @brief Returns the x column of this matrix.
   *
   * @return The x column of this matrix as a vec3<T>.
   */
  vec3<T> get_xcol() const {
    return vec3<T>(xx, yx, zx);
  }

  /**
   * @brief Returns the y column of this matrix.
   *
   * @return The y column of this matrix as a vec3<T>.
   */
  vec3<T> get_ycol() const {
    return vec3<T>(xy, yy, zy);
  }

  /**
   * @brief Returns the z column of this matrix.
   *
   * @return The z column of this matrix as a vec3<T>.
   */
  vec3<T> get_zcol() const {
    return vec3<T>(xz, yz, zz);
  }

  /**
   * @brief Returns the i-th column of this matrix.
   *
   * @param i The column to return, with 0 being the x column, 1 the y column and 2 the z column.
   * @return The i-th column of this matrix as a vec3<T>.
   */
  vec3<T> get_col(int i) const {
    return vec3<T>(*(&xx + i), *(&xx + 3 + i), *(&xx + 6 + i));
  }

  /**
   * @brief Returns the i-th column of this matrix.
   *
   * @param i The column to return, with 0 being the x column, 1 the y column and 2 the z column.
   * @return The i-th column of this matrix as a vec3<T>.
   */
  vec3<T> get_col(size_t i) const {
    return vec3<T>(*(&xx + i), *(&xx + 3 + i), *(&xx + 6 + i));
  }

  /**
   * @brief Computes the inverse of this matrix.
   *
   * If the determinant of this matrix is effectively zero, the inverse is not
   * computed and the zero matrix is returned instead. This is a choice. Why not!
   *
   * @return The inverse of this matrix.
   */
  mat9<T> get_inverse() const {
    double det = xx * (yy * zz - zy * yz) - xy * (yx * zz - yz * zx) + xz * (yx * zy - yy * zx);
    double invdet;
    if (fabs(det) < 1e-20) {
      invdet = 0.0;
    } // this is a choice. Why not!
    else {
      invdet = 1.0 / det;
    }
    return mat9<T>((yy * zz - zy * yz) * invdet, -(xy * zz - xz * zy) * invdet, (xy * yz - xz * yy) * invdet,
                   -(yx * zz - yz * zx) * invdet, (xx * zz - xz * zx) * invdet, -(xx * yz - yx * xz) * invdet,
                   (yx * zy - zx * yy) * invdet, -(xx * zy - zx * xy) * invdet, (xx * yy - yx * xy) * invdet);
  }

  /**
   * @brief Computes the eigenvectors and eigenvalues of this symmetric matrix.
   *
   * This function uses the Jacobi eigenvalue algorithm to compute the
   * eigenvectors and eigenvalues of this symmetric matrix. The eigenvectors are
   * stored as columns in the matrix V and the eigenvalues are stored in the
   * vector D.
   *
   * @param V The matrix in which to store the eigenvectors.
   * @param D The vector in which to store the eigenvalues.
   * @return The number of rotations performed by the algorithm or -1 if the
   * algorithm did not converge.
   */
  int sym_eigen(mat9<T> &V, vec3<T> &D) const {

    /// Compute eigenvectors (stored as columns in V) and corresponding eigenvalues (D)
    /// by assuming the matrix is double and symmetric
    /// See section 11.1 of Numerical Recipes in C for more information.

    int rot = 0;
    // double tresh;
    vec3<T> B;
    vec3<T> Z;

    // Save the input matrix in orig, use new matrix inp
    mat9<T> A = *this;
    // Set vectors to the identity matrix
    V.setIdentity();
    // Set B and D values to the diagonal of the input matrix
    for (size_t i = 0; i < 3; i++) { B[i] = D[i] = A[i * 3 + i]; }

    // Rotate until off-diagonal elements of input matrix are zero
    for (int sweep = 0; sweep++ < 50;) {
      double sum = fabs(A[0 * 3 + 1]) + fabs(A[0 * 3 + 2]) + fabs(A[1 * 3 + 2]);
      double thresh;

      if (fabs(sum) < 1.0e-15) { return rot; }

      thresh = (sweep < 4) ? sum * 0.2 / 9.0 : 0.0; // First three sweeps?

      for (int p = 0; p < 2; p++)
        for (int q = p + 1; q < 3; q++) {
          double g = 100.0 * fabs(A[p * 3 + q]);

          // After 4 sweeps, skip the rotation if the
          // off-diagonal element is small.
          if ((sweep > 4) && (g < 1.0e-15)) A[p * 3 + q] = 0.0;
          else if (fabs(A[p * 3 + q]) > thresh) {
            double h = D[q] - D[p];
            double c, s, t; // cosine, sine, tangent of rotation angle
            double tau;

            if (g < 1.0e-20) t = A[p * 3 + q] / h;
            else {
              double theta = 0.5 * h / A[p * 3 + q];
              t            = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
              if (theta < 0.0) t = -t;
            }

            c   = 1.0 / sqrt(1.0 + t * t); // cosine of rotation angle
            s   = t * c;                   // sine of rotation angle
            tau = s / (1.0 + c);

            h = t * A[p * 3 + q];
            Z[p] -= h;
            Z[q] += h;
            D[p] -= h;
            D[q] += h;
            A[p * 3 + q] = 0.0;

            // case of rotations 0 <= j < p-1
            for (int j = 0; j <= p - 1; j++) {
              g            = A[j * 3 + p];
              h            = A[j * 3 + q];
              A[j * 3 + p] = g - s * (h + g * tau);
              A[j * 3 + q] = h + s * (g - h * tau);
            }

            // case of rotations p < j < q
            for (int j = p + 1; j < q; j++) {
              g            = A[p * 3 + j];
              h            = A[j * 3 + q];
              A[p * 3 + j] = g - s * (h - g * tau);
              A[j * 3 + q] = h + s * (g - h * tau);
            }

            // case of rotations q < j < 3
            for (int j = q + 1; j < 3; j++) {
              g            = A[p * 3 + j];
              h            = A[q * 3 + j];
              A[p * 3 + j] = g - s * (h + g * tau);
              A[q * 3 + j] = h + s * (g - h * tau);
            }

            // Set the eigen vectors
            for (int j = 0; j < 3; j++) {
              g            = V[j * 3 + p];
              h            = V[j * 3 + q];
              V[j * 3 + p] = g - s * (h + g * tau);
              V[j * 3 + q] = h + s * (g - h * tau);
            }
            rot++;
          }
        }

      // Set the eigen values
      B += Z;
      D = B;
      Z.set(0.0, 0.0, 0.0);
    }
    return -1; // Non-normal return - too many rotations
  }

  /**
   * @brief Computes the eigenvectors and eigenvalues of this symmetric matrix, sorted in descending order.
   *
   * This function calls sym_eigen() to compute the eigenvectors and eigenvalues of this symmetric matrix and then
   * sorts them in descending order.
   *
   * @param V The matrix in which to store the eigenvectors.
   * @param D The vector in which to store the eigenvalues.
   */
  void sorted_sym_eigen(mat9<T> &V, vec3<T> &D) {
    this->sym_eigen(V, D);
    // sorting (descending order) bubble sort
    if (D.x < D.y) {
      std::swap(D.x, D.y);
      std::swap(V.xx, V.xy);
      std::swap(V.yx, V.yy);
      std::swap(V.zx, V.zy);
    }
    if (D.y < D.z) {
      std::swap(D.y, D.z);
      std::swap(V.xy, V.xz);
      std::swap(V.yy, V.yz);
      std::swap(V.zy, V.zz);
    }
    if (D.x < D.y) {
      std::swap(D.x, D.y);
      std::swap(V.xx, V.xy);
      std::swap(V.yx, V.yy);
      std::swap(V.zx, V.zy);
    }
  }

  /**
   * @brief Compares this matrix with another matrix for equality.
   *
   * @param other The matrix to compare with this matrix.
   * @return true if all elements of the matrices are equal, false otherwise.
   */
  bool operator==(const mat9<T> &other) const {
    return (this->xx == other.xx && this->xy == other.xy && this->xz == other.xz && this->yx == other.yx &&
            this->yy == other.yy && this->yz == other.yz && this->zx == other.zx && this->zy == other.zy &&
            this->zz == other.zz);
  }

  /**
   * @brief Compares this matrix with another matrix for inequality.
   *
   * @param other The matrix to compare with this matrix.
   * @return true if any elements of the matrices are not equal, false otherwise.
   */
  bool operator!=(const mat9<T> &other) const {
    return !(*this == other);
  }

  friend std::ostream &operator<<(std::ostream &pStr, const mat9 &M) {
    return (pStr << M.xx << ' ' << M.xy << ' ' << M.xz << ' ' << M.yx << ' ' << M.yy << ' ' << M.yz << ' ' << M.zx
                 << ' ' << M.zy << ' ' << M.zz);
  }

  friend std::istream &operator>>(std::istream &pStr, mat9 &M) {
    return (pStr >> M.xx >> M.xy >> M.xz >> M.yx >> M.yy >> M.yz >> M.zx >> M.zy >> M.zz);
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
    auto fmt = [&](T v) {
      std::ostringstream oss;
      oss << std::setprecision(precision);
      if (scientific) oss << std::scientific;
      else oss << std::fixed;
      oss << v;
      return oss.str();
    };

    std::string s[3][3] = {{fmt(xx), fmt(xy), fmt(xz)}, {fmt(yx), fmt(yy), fmt(yz)}, {fmt(zx), fmt(zy), fmt(zz)}};

    // --- Step 2: compute column widths ---
    size_t width[3] = {0, 0, 0};
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i) width[j] = std::max(width[j], s[i][j].size());

    // --- Step 3: big brackets ---
    std::string topLeft = (compact ? "⎡" : "⎡ "), middleLeft = (compact ? "⎢" : "⎢ "),
                bottomLeft = (compact ? "⎣" : "⎣ ");
    std::string topRight = (compact ? "⎤" : " ⎤"), middleRight = (compact ? "⎥" : " ⎥"),
                bottomRight = (compact ? "⎦" : " ⎦");

    if (coloredBrackets) {
      topLeft     = blue + topLeft + reset;
      middleLeft  = blue + middleLeft + reset;
      bottomLeft  = blue + bottomLeft + reset;
      topRight    = blue + topRight + reset;
      middleRight = blue + middleRight + reset;
      bottomRight = blue + bottomRight + reset;
    }

    // --- Step 4: print each row ---
    for (int i = 0; i < 3; ++i) {
      std::ostringstream row;
      for (int j = 0; j < 3; ++j) {
        row << std::setw(static_cast<int>(width[j])) << s[i][j];
        if (j < 2) row << sep;
      }

      std::string lbr = (i == 0) ? topLeft : (i == 2 ? bottomLeft : middleLeft);
      std::string rbr = (i == 0) ? topRight : (i == 2 ? bottomRight : middleRight);

      std::cout << lbr << row.str() << rbr << "\n";
    }
  }
};

// predefined typedefs
typedef mat9<double> mat9r;
typedef mat9<float> mat9f;
typedef mat9<int> mat9i;
typedef mat9<unsigned int> mat9ui;
typedef mat9<bool> mat9b;

/**
 * @brief Computes the dyadic product (tensorial product or otimes) of two vectors.
 *
 * @param a The first vector.
 * @param b The second vector.
 * @return The dyadic product of the two vectors.
 *
 * @note This function returns a matrix in which the element at row i and
 *       column j is the product of the i-th element of the first vector
 *       and the j-th element of the second vector.
 */
template <typename U> mat9<U> dyadic_product(const vec3<U> &a, const vec3<U> &b) {
  return mat9<U>(a.x * b.x, a.x * b.y, a.x * b.z, a.y * b.x, a.y * b.y, a.y * b.z, a.z * b.x, a.z * b.y, a.z * b.z);
}

/**
 * @brief Computes the covariance matrix from a set of points.
 *
 * @param points The set of points.
 * @return The covariance matrix.
 *
 * @note The points are assumed to have the same type as the matrix
 *       elements.
 */
template <typename T> mat9<T> CovarianceMatrix(std::vector<vec3<T>> &points) {
  vec3<T> mu;
  mat9<T> C;

  // loop over the points to find the mean point location
  size_t nbP = 0;
  for (size_t p = 0; p < points.size(); p++) {
    mu += points[p];
    nbP++;
  }

  if (nbP == 0) {
    std::cerr << "@getCovarianceMatrix, no points!\n";
    return C;
  }
  mu /= (T)nbP;

  // loop over the points again to build the covariance matrix.  Note that we only have
  // to build terms for the upper trianglular portion since the matrix is symmetric
  T cxx = 0.0, cxy = 0.0, cxz = 0.0, cyy = 0.0, cyz = 0.0, czz = 0.0;

  for (size_t p = 0; p < points.size(); p++) {
    cxx += points[p].x * points[p].x - mu.x * mu.x;
    cxy += points[p].x * points[p].y - mu.x * mu.y;
    cxz += points[p].x * points[p].z - mu.x * mu.z;
    cyy += points[p].y * points[p].y - mu.y * mu.y;
    cyz += points[p].y * points[p].z - mu.y * mu.z;
    czz += points[p].z * points[p].z - mu.z * mu.z;
  }

  // now build the covariance matrix
  C.xx = cxx;
  C.xy = cxy;
  C.xz = cxz;
  C.yx = cxy;
  C.yy = cyy;
  C.yz = cyz;
  C.zx = cxz;
  C.zy = cyz;
  C.zz = czz;

  return C;
}

#endif /* end of include guard: MAT9_HPP */

#if 0
#include <iostream>

int main (int argc, char const *argv[])
{
  mat9i M(1,2,3,4,5,6,7,8,9);
  //M[4] = 0;
  M.at(0,0) = 12;
  std::cout << M[0] << '\n';
  std::cout << M << '\n';
  
  return 0;
}

#endif

#include "doctest.h"
#include "mat9.hpp"

TEST_CASE("mat9") {
  mat9r a, b;
  b.set_diag(2.0, 2.0, 2.0);
  a.setIdentity();
  REQUIRE(a == mat9r::unit());
  REQUIRE(a + a == b);
  REQUIRE(a * a == a);
}

TEST_CASE("CovarianceMatrix") {
  std::vector<vec3r> points = {
    {1.0, 2.0, 3.0},
    {2.0, 3.0, 4.0},
    {3.0, 4.0, 5.0}
  };

  mat9r expected;
  expected.xx = 2.0;
  expected.xy = 2.0;
  expected.xz = 2.0;
  expected.yx = 2.0;
  expected.yy = 2.0;
  expected.yz = 2.0;
  expected.zx = 2.0;
  expected.zy = 2.0;
  expected.zz = 2.0;

  mat9r result = CovarianceMatrix(points);

  REQUIRE(result == expected);
}

// mat9 is a 3x3 matrix; the 9-arg constructor is row-major (xx,xy,xz, yx,yy,yz, zx,zy,zz).
TEST_CASE("mat9 construction and setters") {
  CHECK(mat9r::zero() == mat9r(0, 0, 0, 0, 0, 0, 0, 0, 0));
  CHECK(mat9r::one() == mat9r(1, 1, 1, 1, 1, 1, 1, 1, 1));
  CHECK(mat9r::unit() == mat9r(1, 0, 0, 0, 1, 0, 0, 0, 1));

  // set(col1, col2, col3) fills the matrix by columns
  mat9r byCols;
  byCols.set(vec3r(1, 4, 7), vec3r(2, 5, 8), vec3r(3, 6, 9));
  CHECK(byCols == mat9r(1, 2, 3, 4, 5, 6, 7, 8, 9));

  // set_col replaces a single column
  mat9r m = mat9r::zero();
  m.set_col(0, vec3r(1, 4, 7));
  CHECK(m.xx == 1.0);
  CHECK(m.yx == 4.0);
  CHECK(m.zx == 7.0);

  mat9r z(5);
  z.setZero();
  CHECK(z == mat9r::zero());
}

TEST_CASE("mat9 arithmetic") {
  mat9r d;
  d.set(2, 0, 0, 0, 3, 0, 0, 0, 4);

  CHECK(d + d == mat9r(4, 0, 0, 0, 6, 0, 0, 0, 8));
  CHECK(d - d == mat9r::zero());
  CHECK(-d == mat9r(-2, 0, 0, 0, -3, 0, 0, 0, -4));
  CHECK(d * 2.0 == mat9r(4, 0, 0, 0, 6, 0, 0, 0, 8));
  CHECK(2.0 * d == mat9r(4, 0, 0, 0, 6, 0, 0, 0, 8));
  CHECK(d / 2.0 == mat9r(1, 0, 0, 0, 1.5, 0, 0, 0, 2));

  mat9r acc = d;
  acc += d;
  CHECK(acc == mat9r(4, 0, 0, 0, 6, 0, 0, 0, 8));
  acc -= d;
  CHECK(acc == d);
  acc *= 2.0;
  CHECK(acc == mat9r(4, 0, 0, 0, 6, 0, 0, 0, 8));
  acc /= 2.0;
  CHECK(acc == d);
}

TEST_CASE("mat9 det, trace, transpose") {
  mat9r d;
  d.set(2, 0, 0, 0, 3, 0, 0, 0, 4);
  CHECK(d.det() == doctest::Approx(24.0));
  CHECK(d.trace() == doctest::Approx(9.0));

  mat9r b(1, 2, 3, 4, 5, 6, 7, 8, 9);
  CHECK(b.transposed() == mat9r(1, 4, 7, 2, 5, 8, 3, 6, 9));
  mat9r t = b;
  t.transpose();
  CHECK(t == mat9r(1, 4, 7, 2, 5, 8, 3, 6, 9));
}

TEST_CASE("mat9 products and inverse") {
  mat9r I = mat9r::unit();
  mat9r d;
  d.set(2, 0, 0, 0, 3, 0, 0, 0, 4);

  CHECK(I * d == d);
  CHECK(d * d == mat9r(4, 0, 0, 0, 9, 0, 0, 0, 16));

  // matrix-vector product
  CHECK(d * vec3r(1, 1, 1) == vec3r(2, 3, 4));
  CHECK(I * vec3r(5, 6, 7) == vec3r(5, 6, 7));

  // A * A^-1 == I (diagonal inverse has an inexact 1/3 entry => use Approx)
  mat9r prod = d * d.get_inverse();
  CHECK(prod.xx == doctest::Approx(1.0));
  CHECK(prod.yy == doctest::Approx(1.0));
  CHECK(prod.zz == doctest::Approx(1.0));
  CHECK(prod.xy == doctest::Approx(0.0));
  CHECK(prod.zx == doctest::Approx(0.0));
}

TEST_CASE("mat9 tensor helpers") {
  mat9r d;
  d.set(2, 0, 0, 0, 3, 0, 0, 0, 4);

  CHECK(inner_product(mat9r::unit(), mat9r::unit()) == doctest::Approx(3.0));
  CHECK(hadamard_product(d, mat9r::one()) == d);

  // spheric + deviatoric == original
  CHECK(spheric(d) == mat9r(3, 0, 0, 0, 3, 0, 0, 0, 3));
  CHECK(deviatoric(d) == mat9r(-1, 0, 0, 0, 0, 0, 0, 0, 1));
  CHECK(spheric(d) + deviatoric(d) == d);
}

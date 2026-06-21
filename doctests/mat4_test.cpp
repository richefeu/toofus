#include "doctest.h"
#include "mat4.hpp"

TEST_CASE("mat4") {
  mat4r a, b;
  b.set_diag(2.0, 2.0);
  a = mat4r::unit();
  REQUIRE( a == mat4r::unit() );
  REQUIRE( a + a == b );
  REQUIRE( a * a == a );
}

// mat4 is a 2x2 matrix (xx, xy, yx, yy) stored row-major.
TEST_CASE("mat4 construction and special matrices") {
  CHECK(mat4r::unit() == mat4r(1, 0, 0, 1));
  CHECK(mat4r::zero() == mat4r(0, 0, 0, 0));

  mat4r d;
  d.set_diag(2.0, 3.0);
  CHECK(d == mat4r(2, 0, 0, 3));

  mat4r id;
  id.setIdentity();
  CHECK(id == mat4r::unit());

  mat4r z(5, 5, 5, 5);
  z.setZero();
  CHECK(z == mat4r::zero());
}

TEST_CASE("mat4 arithmetic") {
  mat4r a(1, 2, 3, 4), b(5, 6, 7, 8);
  CHECK(a + b == mat4r(6, 8, 10, 12));
  CHECK(b - a == mat4r(4, 4, 4, 4));
  CHECK(-a == mat4r(-1, -2, -3, -4));
  CHECK(a * 2.0 == mat4r(2, 4, 6, 8));
  CHECK(2.0 * a == mat4r(2, 4, 6, 8));
  CHECK(a / 2.0 == mat4r(0.5, 1.0, 1.5, 2.0));
}

TEST_CASE("mat4 det, trace, transpose, normSup") {
  mat4r a(1, 2, 3, 4);
  CHECK(a.det() == doctest::Approx(-2.0));
  CHECK(a.trace() == doctest::Approx(5.0));
  CHECK(a.normSup() == doctest::Approx(4.0));
  CHECK(a.transposed() == mat4r(1, 3, 2, 4));

  mat4r t(1, 2, 3, 4);
  t.transpose();
  CHECK(t == mat4r(1, 3, 2, 4));
}

TEST_CASE("mat4 products and inverse") {
  mat4r I = mat4r::unit();
  mat4r a(1, 2, 3, 4);

  CHECK(I * a == a);
  CHECK(a * I == a);

  // matrix-vector product
  CHECK(a * vec2r(1, 1) == vec2r(3, 7));
  CHECK(I * vec2r(2, 5) == vec2r(2, 5));

  // A * A^-1 == I
  mat4r prod = a * a.get_inverse();
  CHECK(prod.xx == doctest::Approx(1.0));
  CHECK(prod.xy == doctest::Approx(0.0));
  CHECK(prod.yx == doctest::Approx(0.0));
  CHECK(prod.yy == doctest::Approx(1.0));

  // in-place inverse returns true and matches get_inverse
  mat4r b(1, 2, 3, 4);
  CHECK(b.inverse() == true);
  CHECK(b == a.get_inverse());

  // a singular matrix cannot be inverted in place
  mat4r singular(1, 1, 1, 1);
  CHECK(singular.inverse() == false);
}



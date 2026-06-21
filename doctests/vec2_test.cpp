#include "doctest.h"
#include "vec2.hpp"

TEST_CASE("vec2") {
  vec2r V;
  REQUIRE(V.x == 0.0);
  REQUIRE(norm(V) == 0.0);
  V.x = 1234.5;
  V.normalize();
  REQUIRE(norm(V) == 1.0);
  REQUIRE(cross(V, V) == 0.0);
  vec2r a(1., 0.);
  vec2r b(0., 1.);
  REQUIRE(a == a);
  REQUIRE(a != b);
  REQUIRE(cross(a, b) == 1.0);
  REQUIRE(cross(a, b) == 1.0);
}

TEST_CASE("vec2 constructors, accessors and setters") {
  vec2r a(2.0, 3.0);
  CHECK(a.x == 2.0);
  CHECK(a.y == 3.0);
  CHECK(a[0] == 2.0);
  CHECK(a[1] == 3.0);

  vec2r copy(a);
  CHECK(copy == a);

  CHECK(vec2r::unit_x() == vec2r(1, 0));
  CHECK(vec2r::unit_y() == vec2r(0, 1));
  CHECK(vec2r::one() == vec2r(1, 1));
  CHECK(vec2r::zero() == vec2r(0, 0));

  vec2r s;
  s.set(4.0, 5.0);
  CHECK(s == vec2r(4, 5));
  s.set(7.0);
  CHECK(s == vec2r(7, 7));
  s.reset();
  CHECK(s == vec2r(0, 0));
  CHECK(s.isnull());
}

TEST_CASE("vec2 arithmetic") {
  vec2r a(1, 2), b(3, 4);
  CHECK(a + b == vec2r(4, 6));
  CHECK(b - a == vec2r(2, 2));
  CHECK(-a == vec2r(-1, -2));
  CHECK(a * 2.0 == vec2r(2, 4));
  CHECK(2.0 * a == vec2r(2, 4));
  CHECK(b / 2.0 == vec2r(1.5, 2.0));

  vec2r c = a;
  c += b;
  CHECK(c == vec2r(4, 6));
  c -= b;
  CHECK(c == a);
  c *= 3.0;
  CHECK(c == vec2r(3, 6));
  c /= 3.0;
  CHECK(c == a);
}

TEST_CASE("vec2 products, norms and helpers") {
  vec2r a(1, 0), b(0, 1);
  CHECK((a * b) == 0.0); // dot product
  CHECK((a * a) == 1.0);
  CHECK(cross(a, b) == 1.0);
  CHECK(determinant(a, b) == 1.0);

  CHECK(component_product(vec2r(2, 3), vec2r(4, 5)) == vec2r(8, 15));
  CHECK(component_min(vec2r(2, 5), vec2r(4, 1)) == vec2r(2, 1));
  CHECK(component_max(vec2r(2, 5), vec2r(4, 1)) == vec2r(4, 5));
  CHECK(component_abs(vec2r(-2, 3)) == vec2r(2, 3));

  vec2r v(3, 4);
  CHECK(norm2(v) == 25.0);
  CHECK(norm(v) == doctest::Approx(5.0));
  CHECK(v.length() == doctest::Approx(5.0));
  CHECK(v.normSup() == 4.0);

  vec2r n = v.normalized(); // const: leaves v untouched
  CHECK(norm(n) == doctest::Approx(1.0));
  CHECK(v == vec2r(3, 4));

  vec2r w(3, 4);
  double len = w.normalize(); // returns the length before normalization
  CHECK(len == doctest::Approx(5.0));
  CHECK(norm(w) == doctest::Approx(1.0));
}

TEST_CASE("vec2 lerp, quarter turns and ordering") {
  CHECK(lerp(0.5, vec2r(0, 0), vec2r(2, 4)) == vec2r(1, 2));
  CHECK(lerp(0.0, vec2r(1, 1), vec2r(9, 9)) == vec2r(1, 1));
  CHECK(lerp(1.0, vec2r(1, 1), vec2r(9, 9)) == vec2r(9, 9));

  vec2r r(2, 3);
  CHECK(r.quarterRightTurned() == vec2r(3, -2));
  CHECK(r.quarterLeftTurned() == vec2r(-3, 2));

  vec2r rr(2, 3);
  rr.quarterRightTurn();
  CHECK(rr == vec2r(3, -2));
  vec2r rl(2, 3);
  rl.quarterLeftTurn();
  CHECK(rl == vec2r(-3, 2));

  // lexicographic operator< (by x, then y)
  CHECK(vec2r(1, 2) < vec2r(2, 0));
  CHECK(vec2r(1, 2) < vec2r(1, 3));
  CHECK_FALSE(vec2r(1, 2) < vec2r(1, 2));
  CHECK_FALSE(vec2r(2, 0) < vec2r(1, 9));
}

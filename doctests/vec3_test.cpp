#include "doctest.h"
#include "vec3.hpp"

TEST_CASE("vec3") {
  vec3r V;
  REQUIRE(V.x == 0.0);
  REQUIRE(norm(V) == 0.0);
  V.x = 1234.5;
  V.normalize();
  REQUIRE(norm(V) == 1.0);
  REQUIRE(norm(cross(V, V)) == 0.0);
  vec3r a(1., 0., 0.);
  vec3r b(0., 1., 0.);
  REQUIRE(a == a);
  REQUIRE(a != b);
  REQUIRE(cross(a, b).z == 1.0);
  REQUIRE(norm(cross(a, b)) == 1.0);
}

TEST_CASE("vec3 constructors, accessors and setters") {
  vec3r a(2.0, 3.0, 4.0);
  CHECK(a[0] == 2.0);
  CHECK(a[1] == 3.0);
  CHECK(a[2] == 4.0);

  vec3r copy(a);
  CHECK(copy == a);

  CHECK(vec3r::unit_x() == vec3r(1, 0, 0));
  CHECK(vec3r::unit_y() == vec3r(0, 1, 0));
  CHECK(vec3r::unit_z() == vec3r(0, 0, 1));
  CHECK(vec3r::one() == vec3r(1, 1, 1));
  CHECK(vec3r::zero() == vec3r(0, 0, 0));

  vec3r s;
  s.set(1.0, 2.0, 3.0);
  CHECK(s == vec3r(1, 2, 3));
  s.set(5.0);
  CHECK(s == vec3r(5, 5, 5));
  s.reset();
  CHECK(s == vec3r(0, 0, 0));
  CHECK(s.isnull());
}

TEST_CASE("vec3 arithmetic") {
  vec3r a(1, 2, 3), b(4, 5, 6);
  CHECK(a + b == vec3r(5, 7, 9));
  CHECK(b - a == vec3r(3, 3, 3));
  CHECK(-a == vec3r(-1, -2, -3));
  CHECK(a * 2.0 == vec3r(2, 4, 6));
  CHECK(2.0 * a == vec3r(2, 4, 6));
  CHECK(b / 2.0 == vec3r(2, 2.5, 3));

  vec3r c = a;
  c += b;
  CHECK(c == vec3r(5, 7, 9));
  c -= b;
  CHECK(c == a);
  c *= 2.0;
  CHECK(c == vec3r(2, 4, 6));
  c /= 2.0;
  CHECK(c == a);
}

TEST_CASE("vec3 products and norms") {
  vec3r x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
  CHECK((x * y) == 0.0); // dot
  CHECK((x * x) == 1.0);
  CHECK(dot(x, y) == 0.0);

  CHECK(cross(x, y) == z);
  CHECK(cross(y, z) == x);
  CHECK(cross(z, x) == y);
  CHECK((x ^ y) == z); // operator^ is the cross product

  CHECK(component_product(vec3r(2, 3, 4), vec3r(1, 2, 3)) == vec3r(2, 6, 12));
  CHECK(component_min(vec3r(2, 5, 1), vec3r(4, 1, 3)) == vec3r(2, 1, 1));
  CHECK(component_max(vec3r(2, 5, 1), vec3r(4, 1, 3)) == vec3r(4, 5, 3));
  CHECK(component_abs(vec3r(-2, 3, -4)) == vec3r(2, 3, 4));

  vec3r v(0, 3, 4);
  CHECK(norm2(v) == 25.0);
  CHECK(norm(v) == doctest::Approx(5.0));
  CHECK(v.length() == doctest::Approx(5.0));
  CHECK(v.normSup() == 4.0);

  CHECK(lerp(0.5, vec3r(0, 0, 0), vec3r(2, 4, 6)) == vec3r(1, 2, 3));
}

TEST_CASE("vec3 normalize variants all yield a unit vector") {
  {
    vec3r v(0, 3, 4);
    double len = v.normalize();
    CHECK(len == doctest::Approx(5.0));
    CHECK(norm(v) == doctest::Approx(1.0));
  }
  {
    vec3r v(1, 2, 2);
    double len = v.normalizeTested();
    CHECK(len == doctest::Approx(3.0));
    CHECK(norm(v) == doctest::Approx(1.0));
  }
  {
    vec3r v(1, 2, 2);
    double len = v.normalizeQuotientAlgo();
    CHECK(len == doctest::Approx(3.0));
    CHECK(norm(v) == doctest::Approx(1.0));
  }
  {
    // normalized() returns a unit copy without mutating the source
    vec3r v(0, 3, 4);
    vec3r n = v.normalized();
    CHECK(norm(n) == doctest::Approx(1.0));
    CHECK(v == vec3r(0, 3, 4));
  }
}

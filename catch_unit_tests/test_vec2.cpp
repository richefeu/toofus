#include "catch.hpp"
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

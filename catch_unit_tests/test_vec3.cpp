#include "catch.hpp"
#include "vec3.hpp"

TEST_CASE("vec3") {
  vec3r V;
  REQUIRE(V.x == 0.0);
  REQUIRE(norm(V) == 0.0);
  V.x = 1234.;
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

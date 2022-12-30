#include "AABB.hpp"
#include "catch.hpp"

TEST_CASE("AABB") {
  AABB aabb1(vec3r(1.0, 1.0, 1.0), vec3r(2.0, 2.0, 2.0));
  AABB aabb2(vec3r(3.0, 3.0, 3.0), vec3r(4.0, 4.0, 4.0));
  REQUIRE(aabb1.intersect(aabb2) == false);
  aabb1.enlarge(1.1);
  REQUIRE(aabb1.intersect(aabb2) == true);
}

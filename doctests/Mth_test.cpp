#include "doctest.h"
#include "Mth.hpp"

TEST_CASE("Mth") {
  REQUIRE(Mth::deg2rad * Mth::rad2deg == 1.0);
  REQUIRE(Mth::sign(-1234.87) == -1.0);
  REQUIRE(Mth::sign(1234.87) == 1.0);
  REQUIRE(Mth::DiamondAngle(1.0, 0.0) == 0.0);
  REQUIRE(Mth::DiamondAngle(0.0, 1.0) == 1.0);
  REQUIRE(Mth::DiamondAngle(-1.0, 0.0) == 2.0);
  REQUIRE(Mth::DiamondAngle(0.0, -1.0) == 3.0);
  REQUIRE(Mth::map(10.0, 0.0, 100.0, 20.0, 30.0) == 21.0);
  REQUIRE(Mth::map(-10.0, 0.0, 100.0, 20.0, 30.0) == 19.0);
  REQUIRE(Mth::lerp(-10.0, 10.0, 0.5) == 0.0);
  REQUIRE(Mth::lerp(0.0, 10.0, 0.5) == 5.0);
  REQUIRE(Mth::lerp(10.0, 0.0, 0.5) == 5.0);
  REQUIRE(Mth::lerp(0.0, 10.0, 1.5) == 15.0);
  REQUIRE(Mth::norm(5.0, 0.0, 10.) == 0.5);
  REQUIRE(Mth::norm(-5.0, 0.0, 10.) == -0.5);
}

#include "catch.hpp"
#include "geoPack3D.hpp"

TEST_CASE("geoPack3D") {
  GeoPack3D GP(2.5, 3, 100, 0, 100, 0, 100, 0, 100);
  GP.seedTime();
  GP.exec();
  
  REQUIRE(GP.sample.empty() == false);
}
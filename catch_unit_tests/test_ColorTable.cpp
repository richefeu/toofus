#include "ColorTable.hpp"
#include "catch.hpp"

TEST_CASE("ColorTable") {
  ColorTable ct;
  ct.setTableID(MATLAB_JET);
  ct.setMinMax(-1.0, 1.0);

  colorRGBA col;

  ct.getRGB(-10.0f, &col);
  REQUIRE(col.r == 0);

  ct.getRGB(0.0f, &col);
  REQUIRE(col.r == 145);

  ct.getRGB(1.0f, &col);
  REQUIRE(col.r == 119);

  ct.getRGB(10.0f, &col);
  REQUIRE(col.r == 119);
}

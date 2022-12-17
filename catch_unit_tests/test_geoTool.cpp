#include "catch.hpp"
#include "geoTool.hpp"

TEST_CASE("geoTool") {
  vec3r v1(0.0, 0.0, 0.0);
  vec3r v2(1.0, 0.0, 0.0);
  vec3r v3(0.0, 1.0, 0.0);
  vec3r v4(0.0, 0.0, 1.0);
  REQUIRE(geoTool::volumeTetrahedron(v1, v2, v3, v4) == Approx(0.166666).epsilon(0.0001));

  mat9r I(0.2, 0.05, 0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.2);
  REQUIRE(geoTool::inertiaDivMassTetrahedron(v1, v2, v3, v4) == I);
}

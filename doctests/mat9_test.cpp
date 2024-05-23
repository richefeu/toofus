#include "doctest.h"
#include "mat9.hpp"

TEST_CASE("mat9") {
  mat9r a, b;
  b.set_diag(2.0, 2.0, 2.0);
  a.setIdentity();
  REQUIRE(a == mat9r::unit());
  REQUIRE(a + a == b);
  REQUIRE(a * a == a);
}

TEST_CASE("CovarianceMatrix") {
  std::vector<vec3r> points = {
    {1.0, 2.0, 3.0},
    {2.0, 3.0, 4.0},
    {3.0, 4.0, 5.0}
  };

  mat9r expected;
  expected.xx = 2.0;
  expected.xy = 2.0;
  expected.xz = 2.0;
  expected.yx = 2.0;
  expected.yy = 2.0;
  expected.yz = 2.0;
  expected.zx = 2.0;
  expected.zy = 2.0;
  expected.zz = 2.0;

  mat9r result = CovarianceMatrix(points);

  REQUIRE(result == expected);
}

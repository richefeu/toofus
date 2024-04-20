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

#include "doctest.h"
#include "quat.hpp"

TEST_CASE("quat") {
  quat q;
  REQUIRE(q.v.x == 0.0);
  REQUIRE(q.s == 1.0);
}
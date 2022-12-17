#include "catch.hpp"
#include "fastSort3.hpp"

TEST_CASE("fastSort3") {
  int a = 1, b = 2, c = 3;
  fastSort3<int>(a, b, c);
  REQUIRE(a == 1);
  REQUIRE(b == 2);
  REQUIRE(c == 3);

  a = 1;
  b = 3;
  c = 0;
  fastSort3<int>(a, b, c);
  REQUIRE(a == 0);
  REQUIRE(b == 1);
  REQUIRE(c == 3);

  a = 2;
  b = 1;
  c = 3;
  fastSort3<int>(a, b, c);
  REQUIRE(a == 1);
  REQUIRE(b == 2);
  REQUIRE(c == 3);

  a = 3;
  b = 1;
  c = 2;
  fastSort3<int>(a, b, c);
  REQUIRE(a == 1);
  REQUIRE(b == 2);
  REQUIRE(c == 3);

  a = 3;
  b = 2;
  c = 1;
  fastSort3<int>(a, b, c);
  REQUIRE(a == 1);
  REQUIRE(b == 2);
  REQUIRE(c == 3);
}

#include "catch.hpp"
#include "mat4.hpp"

TEST_CASE("mat4") {
  mat4r a, b;
  b.set_diag(2.0, 2.0);
  a = mat4r::unit();
  REQUIRE( a == mat4r::unit() );
  REQUIRE( a + a == b );
  REQUIRE( a * a == a );
}



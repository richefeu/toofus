#include "doctest.h"
#include "kwParser.hpp"

TEST_CASE("kwParser") {
  kwParser parser(false);
  
  double value = 0.0;
  parser.kwMap["key"] = __GET__(is, value);
  parser.parseString("key 12.34");
  REQUIRE(value == 12.34);
}


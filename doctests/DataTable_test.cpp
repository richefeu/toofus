#include "doctest.h"
#include "DataTable.hpp"

TEST_CASE("DataTable") {
  DataTable dt;
  
  dt.set("key", 0, 0, 4.56);
  dt.set("key", 0, 1, 1.23);
  REQUIRE(dt.get_ngroup() == 2);
  size_t key_id = dt.get_id("key");
  REQUIRE(dt.get(key_id,0,0) == 4.56);
  REQUIRE(dt.get(key_id,0,1) == 1.23);
  REQUIRE(dt.get(key_id,1,0) == 1.23);  
}
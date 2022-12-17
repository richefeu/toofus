#include "catch.hpp"
#include "message.hpp"

TEST_CASE("message") {
  // std::cout << msg::warn() << "This is a unit test warning\n";
  REQUIRE(msg::HumanReadableSeconds(12345) == "3h 25m 45s");
}

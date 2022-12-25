#include "catch.hpp"
#include "consoleProgressBar.hpp"

TEST_CASE("consoleProgressBar") {
  ConsoleProgressBar cpb(100);
  cpb.setTitle("Title: ");
  cpb.setWidth(50);
  cpb.setProgressChar(124);
  cpb.setVoidChar('-');
  cpb.update(10, std::cout);
  REQUIRE(1 == 1);
}

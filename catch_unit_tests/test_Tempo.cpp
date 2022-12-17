#include "Tempo.hpp"
#include "catch.hpp"

TEST_CASE("Tempo") {
  double V = -1.0;
  double E = -1.0;
  Tempo<double> T;
  T.plug(&V);
  T.plug(&E);
  T.set("Range", 0, 1, 1000, 0.001);
  REQUIRE(T.addresses.size() == 2);

  Tempo<double> T2 = T;
  T2.update(0.1);
  REQUIRE(V == 1000.0);
  REQUIRE(E == 1000.0);

  T2.update(1.1);
  REQUIRE(V == 0.001);
  REQUIRE(E == 0.001);
}

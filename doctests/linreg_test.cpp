#include "doctest.h"
#include "linreg.hpp"

TEST_CASE("linreg") {
  std::vector<double> x(20);
  std::vector<double> y(20);

  for (size_t i = 0; i < 20; i++) {
    x[i] = (double)i * 1.0;
    y[i] = 2.0 * x[i] + 3.0;
  }

  linreg *reg = linreg::get();
  reg->run(x, y);
  REQUIRE(reg->orig == 3.0);
  REQUIRE(reg->slope == 2.0);
  REQUIRE(reg->err == 1.0);

  for (size_t i = 0; i < 20; i++) {
    y[i] = 3.0;
  }
  reg->run(x, y);
  REQUIRE(reg->slope == 0.0);
  REQUIRE(reg->orig == 3.0);
  REQUIRE(reg->err == 1.0);

  for (size_t i = 0; i < 20; i++) {
    y[i] = 2.0 * x[i] + 3.0 + 0.01 * sin(10.0 * x[i]);
  }
  reg->run(x, y);
  REQUIRE(reg->slope == doctest::Approx(2.0).epsilon(0.02));
  REQUIRE(reg->orig == doctest::Approx(3.0).epsilon(0.02));
  REQUIRE(reg->err == doctest::Approx(1.0).epsilon(0.02));
}

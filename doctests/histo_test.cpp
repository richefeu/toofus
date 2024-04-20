#include "doctest.h"
#include "histo.hpp"

#include <vector>

TEST_CASE("histo") {
  std::vector<double> v;
  for (size_t i = 0; i < 1000; i++)
    v.push_back((double)i / 1000.0);

  histo H = histo::pdfNumBins(v, 4);
  REQUIRE(H.data[0].Width == 0.999 / 4.0);
}

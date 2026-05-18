#include "weighted_voro2D.hpp"

#include <cstdlib>
#include <iostream>
#include <random>

int main() {

  lv::LaguerreBuilder<double> b;
  
  // Répartition uniforme avec Mersenne Twister
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (size_t i = 0; i < 1500; i++) { b.add_site(dist(rng), dist(rng), 0.0005); }

  b.compute();

  lv::DCEL<double> dcel;
  lv::build_dcel<double>(b, dcel);

  // visual debug
  lv::SVGOptions opt;

  opt.width  = 500;
  opt.height = 500;
  opt.margin = 50;

  opt.draw_cells       = true;
  opt.draw_sites       = false;
  opt.draw_site_labels = false;
  opt.draw_radii       = false;
  opt.draw_limits      = true;

  lv::write_svg(b, "big.svg", opt);

  return 0;
}
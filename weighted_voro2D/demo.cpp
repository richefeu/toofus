#include "weighted_voro2D.hpp"

#include <iostream>

int main() {

  lv::LaguerreBuilder<double> b;

  // clang-format off
  b.sites = {
      {0.25, 0.25, 0.25}, 
      {0.75, 0.25, 0.25}, 
      {0.95, 0.75, 0.25}, 
      {0.25, 0.75, 0.25}, 
      {0.5, 0.5, 0.1},
  };
  // clang-format on

  // b.set_bbox(0.0, 0.0, 1.0, 1.0);
  // b.fit_bbox();
  b.compute_fast2();

  // dcel
  lv::DCEL<double> dcel;
  lv::build_dcel<double>(b, dcel);

  // visual debug
  lv::SVGOptions opt;

  opt.width  = 500;
  opt.height = 500;
  opt.margin = 50;

  opt.draw_cells       = false;
  opt.draw_sites       = true;
  opt.draw_site_labels = true;
  opt.draw_radii       = true;
  opt.draw_limits      = false;
  lv::write_svg_dcel(dcel, b, "debug.svg", opt);

  opt.draw_cells  = true;
  opt.draw_limits = true;
  lv::write_svg(b, "quick.svg", opt);

  return 0;
}
#include "weighted_voro2D.hpp"

#include <fstream>
#include <iostream>

int main() {

  lv::LaguerreBuilder<double> b;

  std::ifstream f("disks.txt");
  double x, y, r;
  while (f.good()) {
    f >> x >> y >> r;
    b.add_site(x, y, r);
  }

  b.add_surrounding_sites(0.21);
  b.compute_fast2();

  lv::write_svg(b, "maple.svg",
                {.width            = 1000,
                 .height           = 1000,
                 .margin           = 30,
                 .draw_cells       = true,
                 .draw_sites       = false,
                 .draw_site_labels = false,
                 .draw_radii       = false,
                 .draw_limits      = true});

  lv::DCEL<double> dcel;
  lv::build_dcel<double>(b, dcel);

  lv::write_svg_dcel(dcel, b, "maple_dcel.svg",
                     {.width            = 2000,
                      .height           = 2000,
                      .margin           = 30,
                      .draw_cells       = false,
                      .draw_sites       = false,
                      .draw_site_labels = false,
                      .draw_radii       = false,
                      .draw_limits      = false,
                      .draw_dcel_labels = false});

  return 0;
}
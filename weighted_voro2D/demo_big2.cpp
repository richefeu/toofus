#include "weighted_voro2D.hpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>

int main() {

  lv::LaguerreBuilder<double> b;

  // Répartition uniforme avec Mersenne Twister
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (size_t i = 0; i < 5000; i++) { b.add_site(dist(rng), dist(rng), 0.0005); }

  // Mesure compute()
  auto t0 = std::chrono::high_resolution_clock::now();
  // b.compute(); // 500:29    2000:259      5000:1493
  // b.compute_fast(); //  500:15    2000:157      5000:926
  b.compute_fast2(); //  500:2    2000:9      5000:19
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "compute()               : " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << " ms\n";

  // Mesure build_true_laguerre_dcel()
  lv::DCEL<double> dcel;
  auto t2 = std::chrono::high_resolution_clock::now();
  lv::build_dcel<double>(b, dcel);
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "build_true_laguerre_dcel: " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
            << " ms\n";

  // Options SVG
  lv::SVGOptions opt;
  opt.width            = 1500;
  opt.height           = 1500;
  opt.margin           = 50;
  opt.draw_cells       = true;
  opt.draw_sites       = false;
  opt.draw_site_labels = false;
  opt.draw_radii       = true;
  opt.draw_limits      = true;

  // Mesure write_svg()
  auto t4 = std::chrono::high_resolution_clock::now();
  lv::write_svg(b, "big.svg", opt);
  auto t5 = std::chrono::high_resolution_clock::now();
  std::cout << "write_svg()             : " << std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count()
            << " ms\n";

  return 0;
}
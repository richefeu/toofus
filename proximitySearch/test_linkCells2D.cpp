#include "linkCells2D.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

struct Grain {
  double xs{0.0}, ys{0.0}, rs{0.0};
  Grain(double x, double y, double r) : xs(x), ys(y), rs(r) {}
};

struct Interaction {
  int i{0}, j{0};
  Interaction(int I, int J) : i(I), j(J) {}
};

int main() {
  std::vector<Grain> grains;

  lc2d::Bounds bounds;
  bounds.xmin = 1e12f;
  bounds.ymin = 1e12f;
  bounds.xmax = -1e12f;
  bounds.ymax = -1e12f;

  std::ifstream f("disks.txt");
  float x, y, r;
  while (f >> x >> y >> r) {
    grains.emplace_back(x, y, r);
    bounds.xmin = std::min(bounds.xmin, x - r);
    bounds.xmax = std::max(bounds.xmax, x + r);
    bounds.ymin = std::min(bounds.ymin, y - r);
    bounds.ymax = std::max(bounds.ymax, y + r);
  }

  double RMAX = 0.02;
  double dstMAX = 0.01;
  //double dstMAX = 0.1; // sur plusieurs horizons

  // O(n^2) ------------------------------------------
  auto t0 = std::chrono::high_resolution_clock::now();
  std::vector<Interaction> inters0;
  for (size_t i = 0; i < grains.size(); i++) {
    for (size_t j = i + 1; j < grains.size(); j++) {
      double dx  = grains[j].xs - grains[i].xs;
      double dy  = grains[j].ys - grains[i].ys;
      double sum = grains[j].rs + grains[i].rs + dstMAX;
      if ((dx * dx + dy * dy) < sum * sum) { inters0.emplace_back(i, j); }
    }
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "brute force : " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << " ms, nbPairs: " << inters0.size() << std::endl;

  // link-cells ---------------------------------------
  auto t2 = std::chrono::high_resolution_clock::now();

  // Les accesseurs sont déclarés une seule fois ici
  auto getX = [&](uint32_t i) { return grains[i].xs; };
  auto getY = [&](uint32_t i) { return grains[i].ys; };
  auto getR = [&](uint32_t i) { return grains[i].rs; };

  lc2d::LinkCells2D grid(bounds, 2 * RMAX + dstMAX, getX, getY, getR);

  grid.build(grains.size());

  std::vector<Interaction> inters;
  grid.forEachNeighborPair(dstMAX, [&](uint32_t i, uint32_t j) { inters.emplace_back(i, j); });

  // étapes supplémentaires nécessaires pour le test de comparaison
  for (auto &interaction : inters) {
    if (interaction.j < interaction.i) { std::swap(interaction.i, interaction.j); }
  }

  // Tri lexicographique : d'abord par i, puis par j
  std::sort(inters.begin(), inters.end(), [](const Interaction &a, const Interaction &b) {
    if (a.i != b.i) {
      return a.i < b.i; // Compare i d'abord
    }
    return a.j < b.j; // Si i est égal, compare j
  });

  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "link cells  : " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
            << " ms, nbPairs: " << inters.size() << std::endl;

  // check ----------------------------------------------
  size_t nbok = 0;
  if (inters.size() == inters0.size()) {
    for (size_t i = 0; i < inters0.size(); i++) {
      if (inters0[i].i == inters[i].i && inters0[i].j == inters[i].j) nbok++;
    }
    std::cout << "nbok = " << nbok << std::endl;
  }

  grid.writeSVG("debug.svg", {.svgWidth = 1000, .drawPairs = true, .dmax = (float)(dstMAX)});

  return 0;
}
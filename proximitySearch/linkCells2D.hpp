#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace lc2d {

struct Vec2 {
  float x = 0.f;
  float y = 0.f;
};

struct Bounds {
  float xmin = 0.f;
  float ymin = 0.f;
  float xmax = 1.f;
  float ymax = 1.f;

  float width() const {
    return xmax - xmin;
  }
  float height() const {
    return ymax - ymin;
  }
};

struct Statistics {
  uint32_t particleCount  = 0;
  uint32_t cellCount      = 0;
  uint32_t maxOccupancy   = 0;
  double averageOccupancy = 0.0;
  uint32_t nonEmptyCells  = 0;
};

// Helper for CTAD
template <class XAcc, class YAcc, class RAcc> class LinkCells2D {
public:
  // ------------------------------------------------------------------ ctor --

  LinkCells2D(const Bounds &b, float cellSize, XAcc getX, YAcc getY, RAcc getR)
      : getX_(std::move(getX)), getY_(std::move(getY)), getR_(std::move(getR)) {
    reset(b, cellSize);
  }

  // ----------------------------------------------------------------- reset --

  void reset(const Bounds &b, float cellSize) {
    bounds_      = b;
    cellSize_    = cellSize;
    invCellSize_ = 1.f / cellSize_;

    nx_ = std::max(1u, uint32_t(std::ceil(bounds_.width() * invCellSize_)));
    ny_ = std::max(1u, uint32_t(std::ceil(bounds_.height() * invCellSize_)));

    cellCount_ = nx_ * ny_;

    counts_.assign(cellCount_, 0);
    offsets_.assign(cellCount_ + 1, 0);
    cursor_.assign(cellCount_, 0);
  }

  // ------------------------------------------------------------ accessors ---

  const Bounds &bounds() const {
    return bounds_;
  }
  uint32_t nx() const {
    return nx_;
  }
  uint32_t ny() const {
    return ny_;
  }
  uint32_t cellCount() const {
    return cellCount_;
  }
  float cellSize() const {
    return cellSize_;
  }

  // ----------------------------------------------------------------- build --

  void build(uint32_t n) {
    particleCount_ = n;

    std::fill(counts_.begin(), counts_.end(), 0);
    particleCell_.resize(n);

    for (uint32_t i = 0; i < n; ++i) {
      const uint32_t c = computeCell(getX_(i), getY_(i));
      particleCell_[i] = c;
      counts_[c]++;
    }

    offsets_[0] = 0;
    for (uint32_t c = 0; c < cellCount_; ++c) offsets_[c + 1] = offsets_[c] + counts_[c];

    ids_.resize(n);
    cursor_ = offsets_;

    for (uint32_t i = 0; i < n; ++i) {
      const uint32_t c   = particleCell_[i];
      ids_[cursor_[c]++] = i;
    }
  }

  // ------------------------------------------------- forEachCandidatePair --

  template <class Callback> void forEachCandidatePair(Callback &&callback) const {
    static constexpr int stencil[5][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};

    for (uint32_t cy = 0; cy < ny_; ++cy) {
      for (uint32_t cx = 0; cx < nx_; ++cx) {
        const uint32_t c0     = cellIndex(cx, cy);
        const uint32_t begin0 = offsets_[c0];
        const uint32_t end0   = offsets_[c0 + 1];

        if (begin0 == end0) continue;

        for (const auto &s : stencil) {
          const int nx = int(cx) + s[0];
          const int ny = int(cy) + s[1];

          if (nx < 0 || ny < 0) continue;
          if (nx >= int(nx_) || ny >= int(ny_)) continue;

          const uint32_t c1     = cellIndex(uint32_t(nx), uint32_t(ny));
          const uint32_t begin1 = offsets_[c1];
          const uint32_t end1   = offsets_[c1 + 1];

          if (begin1 == end1) continue;

          if (c0 == c1) {
            for (uint32_t a = begin0; a < end0; ++a)
              for (uint32_t b = a + 1; b < end0; ++b) callback(ids_[a], ids_[b]);
          } else {
            for (uint32_t a = begin0; a < end0; ++a)
              for (uint32_t b = begin1; b < end1; ++b) callback(ids_[a], ids_[b]);
          }
        }
      }
    }
  }

  // -------------------------------------------------- forEachNeighborPair --

  template <class Callback> void forEachNeighborPair(float dmax, Callback &&callback) const {
    forEachCandidatePair([&](uint32_t i, uint32_t j) {
      const float dx    = getX_(j) - getX_(i);
      const float dy    = getY_(j) - getY_(i);
      const float limit = getR_(i) + getR_(j) + dmax;
      if (dx * dx + dy * dy <= limit * limit) callback(i, j);
    });
  }

  // -------------------------------------------------------------- statistics

  Statistics statistics() const {
    Statistics s;
    s.particleCount = particleCount_;
    s.cellCount     = cellCount_;

    uint64_t sum = 0;
    for (uint32_t c = 0; c < cellCount_; ++c) {
      const uint32_t n = offsets_[c + 1] - offsets_[c];
      if (n > 0) s.nonEmptyCells++;
      s.maxOccupancy = std::max(s.maxOccupancy, n);
      sum += n;
    }

    if (s.nonEmptyCells > 0) s.averageOccupancy = double(sum) / double(s.nonEmptyCells);

    return s;
  }

  // ---------------------------------------------------------------- writeSVG

  struct OptSVG {
    uint32_t svgWidth = 800;
    bool drawGrid     = true;
    bool drawIDs      = false;
    bool drawPairs    = false;
    float dmax        = 0.f;
  };

  void writeSVG(const std::string &filename, const OptSVG &opt = {}) const {
    std::ofstream file(filename);

    const float W  = bounds_.width();
    const float H  = bounds_.height();
    const float sw = W / 1000.0;

    const uint32_t svgH = uint32_t(std::round(opt.svgWidth * H / W));

    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "width=\"" << opt.svgWidth << "\" height=\"" << svgH << "\" "
         << "viewBox=\"" << bounds_.xmin << " " << bounds_.ymin << " " << W << " " << H << "\">\n";

    file << "<style> *{stroke-width: " << 1.5f * W / (float)opt.svgWidth << ";} </style>\n";

    file << "<rect x=\"" << bounds_.xmin << "\" y=\"" << bounds_.ymin << "\" width=\"" << W << "\" height=\"" << H
         << "\" fill=\"white\" stroke=\"black\"/>\n";

    if (opt.drawGrid) {
      file << "<g stroke=\"#cccccc\">\n";

      for (uint32_t ix = 0; ix <= nx_; ++ix) {
        const float x = bounds_.xmin + ix * cellSize_;
        file << "<line x1=\"" << x << "\" y1=\"" << bounds_.ymin << "\" x2=\"" << x << "\" y2=\"" << bounds_.ymax
             << "\"/>\n";
      }
      for (uint32_t iy = 0; iy <= ny_; ++iy) {
        const float y = bounds_.ymin + iy * cellSize_;
        file << "<line x1=\"" << bounds_.xmin << "\" y1=\"" << y << "\" x2=\"" << bounds_.xmax << "\" y2=\"" << y
             << "\"/>\n";
      }
      file << "</g>\n";
    }

    if (opt.drawPairs) {
      file << "<g stroke=\"#66aaff\">\n";
      forEachNeighborPair(opt.dmax, [&](uint32_t i, uint32_t j) {
        file << "<line x1=\"" << getX_(i) << "\" y1=\"" << getY_(i) << "\" x2=\"" << getX_(j) << "\" y2=\"" << getY_(j)
             << "\"/>\n";
      });
      file << "</g>\n";
    }

    file << "<g fill=\"none\" stroke=\"black\">\n";
    for (uint32_t i = 0; i < particleCount_; ++i) {
      file << "<circle cx=\"" << getX_(i) << "\" cy=\"" << getY_(i) << "\" r=\"" << getR_(i) << "\"/>\n";
    }
    file << "</g>\n";

    if (opt.drawIDs) {
      file << "<g font-size=\"0.03\" fill=\"red\">\n";
      for (uint32_t i = 0; i < particleCount_; ++i) {
        file << "<text x=\"" << getX_(i) << "\" y=\"" << getY_(i) << "\">" << i << "</text>\n";
      }
      file << "</g>\n";
    }

    file << "</svg>\n";
  }

private:
  // --------------------------------------------------------------- helpers --

  uint32_t computeCell(float x, float y) const {
    int ix = int((x - bounds_.xmin) * invCellSize_);
    int iy = int((y - bounds_.ymin) * invCellSize_);
    ix     = std::max(0, std::min(ix, int(nx_ - 1)));
    iy     = std::max(0, std::min(iy, int(ny_ - 1)));
    return cellIndex(uint32_t(ix), uint32_t(iy));
  }

  uint32_t cellIndex(uint32_t ix, uint32_t iy) const {
    return ix + iy * nx_;
  }

private:
  // ------------------------------------------------------------ accessors ---
  XAcc getX_;
  YAcc getY_;
  RAcc getR_;

  // ------------------------------------------------------------- geometry ---
  Bounds bounds_;
  float cellSize_     = 1.f;
  float invCellSize_  = 1.f;
  uint32_t nx_        = 1;
  uint32_t ny_        = 1;
  uint32_t cellCount_ = 1;

  // ------------------------------------------------------------ particles ---
  uint32_t particleCount_ = 0;

  std::vector<uint32_t> counts_;
  std::vector<uint32_t> offsets_;
  std::vector<uint32_t> cursor_;
  std::vector<uint32_t> ids_;
  std::vector<uint32_t> particleCell_;
};

// -------------------------------------------------- deduction guide (C++17) --
// Permet d'écrire : lc2d::LinkCells2D grid(bounds, size, getX, getY, getR);
// sans spécifier les types template explicitement.
template <class XAcc, class YAcc, class RAcc>
LinkCells2D(const Bounds &, float, XAcc, YAcc, RAcc) -> LinkCells2D<XAcc, YAcc, RAcc>;

} // namespace lc2d

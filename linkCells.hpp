// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#ifndef LINKCELLS_HPP
#define LINKCELLS_HPP

#include <vector>

#include "AABB.hpp"
#include "vec3.hpp"

enum linkCellOptions {
  PERIODIC_LINKCELLS = 1,
  HALF_CONNECTED_LINKCELLS = 1 << 1,
};

// An axis aligned box used for neighbor tracking
class AABB_Cell {
public:
  std::vector<size_t> bodies;      ///< Holded bodies identifiers
  std::vector<AABB_Cell *> pcells; ///< Surroundind cells (+ current cell)

  AABB_Cell() {}
  ~AABB_Cell() {}
};

class linkCells {
private:
  vec3r factor; ///< Used to convert position to index

public:
  AABB box;                                               ///< Overall surrounding box
  vec3r minSizes;                                         ///< Wanted minimum size (along x, y and z) of the cells
  std::vector<std::vector<std::vector<AABB_Cell>>> cells; ///< Cells that hold only free bodies
  AABB_Cell oversized_bodies;                             ///< A particular cell that hold only 'too big' bodies
  vec3ui N;                                               ///< Number of cells in each direction (x, y, and z)

  /*
  changements :
  1. parametres optionnels pour gerer les cas periodic et demi connectes
  2. parametre optionnel pour ajouter le lien vers oversized_bodies aux pcells
  */

  // Ctor
  linkCells(AABB &Box, vec3r &CellMinSizes, int options = 0) : box(Box), minSizes(CellMinSizes) { init(options); }

  // Dtor
  ~linkCells() {}

  void init(int options = 0) {
    // Partition
    N.x = (unsigned int)floor((box.max.x - box.min.x) / minSizes.x);
    N.y = (unsigned int)floor((box.max.y - box.min.y) / minSizes.y);
    N.z = (unsigned int)floor((box.max.z - box.min.z) / minSizes.z);

    // We've got at least one element per side
    if (N.x < 1)
      N.x = 1;
    if (N.y < 1)
      N.y = 1;
    if (N.z < 1)
      N.z = 1;

    // Cell dimensions
    double dx, dy, dz;
    dx = (box.max.x - box.min.x) / (double)N.x;
    dy = (box.max.y - box.min.y) / (double)N.y;
    dz = (box.max.z - box.min.z) / (double)N.z;

    factor.set((dx > 0.0) ? 1.0 / dx : 0.0, (dy > 0.0) ? 1.0 / dy : 0.0, (dz > 0.0) ? 1.0 / dz : 0.0);

    // Reserve memory
    AABB_Cell A;
    std::vector<AABB_Cell> kvec;
    for (size_t k = 0; k < N.z; ++k) {
      kvec.push_back(A);
    }
    std::vector<std::vector<AABB_Cell>> jvec;
    for (size_t j = 0; j < N.y; ++j)
      jvec.push_back(kvec);
    for (size_t i = 0; i < N.x; ++i)
      cells.push_back(jvec);

    // Link the cells
    size_t ix0, ix1;
    size_t iy0, iy1;
    size_t iz0, iz1;
    for (size_t ix = 0; ix < N.x; ++ix) {
      for (size_t iy = 0; iy < N.y; ++iy) {
        for (size_t iz = 0; iz < N.z; ++iz) {

          if (options & PERIODIC_LINKCELLS) {

            ix0 = (ix > 0) ? ix - 1 : ix;
            ix1 = (ix < N.x - 1) ? ix + 1 : ix;
            iy0 = (iy > 0) ? iy - 1 : iy;
            iy1 = (iy < N.y - 1) ? iy + 1 : iy;
            iz0 = (iz > 0) ? iz - 1 : iz;
            iz1 = (iz < N.z - 1) ? iz + 1 : iz;

          } else {

            ix0 = (ix > 0) ? ix - 1 : N.x - 1;
            ix1 = (ix < N.x - 1) ? ix + 1 : 0;
            iy0 = (iy > 0) ? iy - 1 : N.y - 1;
            iy1 = (iy < N.y - 1) ? iy + 1 : 0;
            iz0 = (iz > 0) ? iz - 1 : N.z - 1;
            iz1 = (iz < N.z - 1) ? iz + 1 : 0;
          }

          /*
if (options & HALF_CONNECTED_LINKCELLS) {

for (size_t iix = ix0; iix <= ix1; ++iix) {
for (size_t iiy = iy0; iiy <= iy1; ++iiy) {
for (size_t iiz = iz0; iiz <= iz1; ++iiz) {
cells[ix][iy][iz].pcells.push_back(&cells[iix][iiy][iiz]);
}
}
}

} else {

for (size_t iix = ix; iix <= ix1; ++iix) {
for (size_t iiy = iy; iiy <= iy1; ++iiy) {
for (size_t iiz = iz; iiz <= iz1; ++iiz) {
cells[ix][iy][iz].pcells.push_back(&cells[iix][iiy][iiz]);
}
}
}

          }
          */

          std::vector<size_t> xxx{ix0, ix, ix1};
          std::vector<size_t> yyy{iy0, iy, iy1};
          std::vector<size_t> zzz{iz0, iz, iz1};

          for (auto iix : xxx) {
            for (auto iiy : yyy) {
              for (auto iiz : zzz) {
                cells[ix][iy][iz].pcells.push_back(&cells[iix][iiy][iiz]);
              }
            }
          }

          cells[ix][iy][iz].pcells.push_back(&oversized_bodies); // any cell 'see' the oversized_bodies
        }
      }
    }
  }

  void clear() {
    for (size_t ix = 0; ix < N.x; ++ix) {
      for (size_t iy = 0; iy < N.y; ++iy) {
        for (size_t iz = 0; iz < N.z; ++iz) {
          cells[ix][iy][iz].bodies.clear();
        }
      }
    }
    oversized_bodies.bodies.clear();
  }

  void add_body(size_t B, vec3r &pos, vec3r diag) {
    int ix, iy, iz;

    if (diag.x > minSizes.x || diag.y > minSizes.y || diag.z > minSizes.z)
      oversized_bodies.bodies.push_back(B);
    else {
      ix = (int)trunc((pos.x - box.min.x) * factor.x);
      iy = (int)trunc((pos.y - box.min.y) * factor.y);
      iz = (int)trunc((pos.z - box.min.z) * factor.z);

      if (ix < 0 || ix >= (int)N.x) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      if (iy < 0 || iy >= (int)N.y) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      if (iz < 0 || iz >= (int)N.z) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      cells[ix][iy][iz].bodies.push_back(B);
    }
  }

  void add_body(size_t B, vec3r &pos, AABB &aabb) {
    int ix, iy, iz;
    vec3r diag = aabb.max - aabb.min;

    if (diag.x > minSizes.x || diag.y > minSizes.y || diag.z > minSizes.z)
      oversized_bodies.bodies.push_back(B);
    else {
      ix = (int)trunc((pos.x - box.min.x) * factor.x);
      iy = (int)trunc((pos.y - box.min.y) * factor.y);
      iz = (int)trunc((pos.z - box.min.z) * factor.z);

      if (ix < 0 || ix >= (int)N.x) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      if (iy < 0 || iy >= (int)N.y) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      if (iz < 0 || iz >= (int)N.z) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      cells[ix][iy][iz].bodies.push_back(B);
    }
  }

  void add_body(size_t B, vec3r &pos) {
    int ix, iy, iz;

    ix = (int)trunc((pos.x - box.min.x) * factor.x);
    iy = (int)trunc((pos.y - box.min.y) * factor.y);
    iz = (int)trunc((pos.z - box.min.z) * factor.z);

    /*
if (ix < 0 || ix >= (int)N.x) {
oversized_bodies.bodies.push_back(B);
return;
}
if (iy < 0 || iy >= (int)N.y) {
oversized_bodies.bodies.push_back(B);
return;
}
if (iz < 0 || iz >= (int)N.z) {
oversized_bodies.bodies.push_back(B);
return;
}
    */
    cells[ix][iy][iz].bodies.push_back(B);
  }
};

#endif /* end of include guard: LINKCELLS_HPP */

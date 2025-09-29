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

#ifndef LINKCELLS_2D_HPP
#define LINKCELLS_2D_HPP

#include <vector>

#include "AABB_2D.hpp"
#include "vec2.hpp"

/// An axis aligned box used for neighbor tracking
class AABB_2D_Cell {
public:
  std::vector<size_t> bodies;         ///< Holded bodies identifiers
  std::vector<AABB_2D_Cell *> pcells; ///< Surroundind cells (+ current cell)

  AABB_2D_Cell() {}
  ~AABB_2D_Cell() {}
};

class linkCells2D {
private:
  vec2r factor; ///< Used to convert position to index

public:
  AABB_2D box;                                  ///< Overall surrounding box
  vec2r minSizes;                               ///< Wanted minimum size (along x, y and z) of the cells
  std::vector<std::vector<AABB_2D_Cell>> cells; ///< Cells that hold only free bodies
  AABB_2D_Cell oversized_bodies;                ///< A particular cell that hold only 'too big' bodies
  vec2ui N;                                     ///< Number of cells in each direction (x, y, and z)

  // Ctor
  linkCells2D(AABB_2D &Box, vec2r &CellMinSizes) : box(Box), minSizes(CellMinSizes) {
    init();
  }

  // Dtor
  ~linkCells2D() {}

  /**
   * Initializes the linkCells2D object by partitioning the overall bounding box
   * into a grid of smaller cells. The function calculates the number of cells
   * in each direction (x and y) based on the minimum cell sizes and ensures
   * there is at least one cell per direction. It then computes the cell
   * dimensions and sets up a conversion factor for mapping positions to cell
   * indices. Memory is allocated for the cell grid, and each cell is linked to
   * its neighboring cells and the oversized bodies cell.
   */
  void init() {
    // Partition
    N.x = (unsigned int)floor((box.max.x - box.min.x) / minSizes.x);
    N.y = (unsigned int)floor((box.max.y - box.min.y) / minSizes.y);

    // We've got at least one element per side
    if (N.x < 1) N.x = 1;
    if (N.y < 1) N.y = 1;

    // Cell dimensions
    double dx, dy;
    dx = (box.max.x - box.min.x) / (double)N.x;
    dy = (box.max.y - box.min.y) / (double)N.y;

    factor.set((dx > 0.0) ? 1.0 / dx : 0.0, (dy > 0.0) ? 1.0 / dy : 0.0);

    // Reserve memory
    AABB_2D_Cell A;
    std::vector<AABB_2D_Cell> jvec;
    for (size_t j = 0; j < N.y; ++j) jvec.push_back(A);
    for (size_t i = 0; i < N.x; ++i) cells.push_back(jvec);

    // Link the cells
    size_t ix0, ix1;
    size_t iy0, iy1;
    for (size_t ix = 0; ix < N.x; ++ix) {
      for (size_t iy = 0; iy < N.y; ++iy) {

        ix0 = (ix > 0) ? ix - 1 : N.x - 1;
        ix1 = (ix < N.x - 1) ? ix + 1 : 0;
        iy0 = (iy > 0) ? iy - 1 : N.y - 1;
        iy1 = (iy < N.y - 1) ? iy + 1 : 0;

        std::vector<size_t> xxx{ix0, ix, ix1};
        std::vector<size_t> yyy{iy0, iy, iy1};

        for (auto iix : xxx) {
          for (auto iiy : yyy) { cells[ix][iy].pcells.push_back(&cells[iix][iiy]); }
        }

        cells[ix][iy].pcells.push_back(&oversized_bodies); // any cell 'see' the oversized_bodies
      }
    }
  }

  /**
   * @brief Clear all cells and oversized bodies.
   *
   * This method clears all cells and oversized bodies of any body. It is
   * typically called after all bodies have been updated.
   */
  void clear() {
    for (size_t ix = 0; ix < N.x; ++ix) {
      for (size_t iy = 0; iy < N.y; ++iy) { cells[ix][iy].bodies.clear(); }
    }
    oversized_bodies.bodies.clear();
  }

  /**
   * @brief Add a body to a cell.
   *
   * This method adds a body to a cell given its position and AABB. If the
   * body is too large for the grid, it is added to the list of oversized
   * bodies.
   *
   * @param B the body ID
   * @param pos the position of the body
   * @param aabb the AABB of the body
   */
  void add_body(size_t B, vec2r &pos, AABB_2D &aabb) {
    int ix, iy;
    // vec2r diag =  aabb.max - aabb.min;
    double surf = fabs((aabb.max.x - aabb.min.x) * (aabb.max.y - aabb.min.y));
    // if (diag.x > minSizes.x || diag.y > minSizes.y)
    if (surf >= fabs(minSizes.x * minSizes.y)) oversized_bodies.bodies.push_back(B);
    else {
      ix = (int)trunc((pos.x - box.min.x) * factor.x);
      iy = (int)trunc((pos.y - box.min.y) * factor.y);

      if (ix < 0 || ix >= (int)N.x) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      if (iy < 0 || iy >= (int)N.y) {
        oversized_bodies.bodies.push_back(B);
        return;
      }
      cells[ix][iy].bodies.push_back(B);
    }
  }

  /**
   * @brief Add a body to a cell based on its position.
   *
   * This method calculates the cell indices for the given position and adds
   * the body to the corresponding cell. It assumes that the position is within
   * the bounds of the overall bounding box.
   *
   * @param B The body ID to add.
   * @param pos The position of the body.
   */
  void add_body(size_t B, vec2r &pos) {
    int ix, iy;

    ix = (int)trunc((pos.x - box.min.x) * factor.x);
    iy = (int)trunc((pos.y - box.min.y) * factor.y);

    cells[ix][iy].bodies.push_back(B);
  }
};

#endif /* end of include guard: LINKCELLS_2D_HPP */

#if 0

int main(int argc, char const *argv[]) { return 0; }

#endif

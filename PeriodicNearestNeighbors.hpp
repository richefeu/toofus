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

#ifndef PERIODIC_NN_HPP
#define PERIODIC_NN_HPP

#include <cmath>
#include <set>
#include <vector>

struct IdPoint {
  size_t id;
  double x, y, z;
  IdPoint(size_t t_id, double t_x, double t_y, double t_z) : id(t_id), x(t_x), y(t_y), z(t_z) {}
};

class PeriodicNearestNeighbors {
public:
  /// Constructs a PeriodicNearestNeighbors object with the given points and gridSize.
  ///
  /// @param t_points The vector of points.
  /// @param t_gridSize The size of the grid.
  ///
  PeriodicNearestNeighbors(const std::vector<IdPoint> &t_points, size_t t_gridSize)
      : points(t_points), gridSize(t_gridSize), numCells(t_gridSize * t_gridSize * t_gridSize),
        gridSizeDouble(static_cast<double>(gridSize)), gridSizeInt(static_cast<int>(gridSize)),
        gridSizeSquared(gridSize * gridSize) {
    buildGrid();
  }

  std::vector<std::vector<size_t>> getNeighbors(double dmax) {
    const size_t numPoints = points.size();
    std::vector<std::vector<size_t>> neighbors(numPoints);

    for (size_t i = 0; i < numPoints; ++i) {
      const IdPoint &point = points[i];
      size_t cellIndex = calculateCellIndex(point);
      std::set<size_t> uniqueNeighborIndices;

      const std::vector<size_t> &neighborIndices = neighborCellIndices[cellIndex];
      for (size_t neighborCellIndex : neighborIndices) {
        for (size_t neighborIndex : grid[neighborCellIndex]) {
          if (neighborIndex <= i)
            continue;

          const IdPoint &neighbor = points[neighborIndex];
          double sqrDistance = calculateSqrDistance(point, neighbor);

          if (sqrDistance <= dmax * dmax) {
            uniqueNeighborIndices.insert(neighborIndex);
          }
        }
      }

      // neighbors[i].reserve(uniqueNeighborIndices.size());
      for (size_t neighborIndex : uniqueNeighborIndices) {
        neighbors[i].push_back(neighborIndex);
      }
    }

    return neighbors;
  }

  std::vector<std::vector<size_t>> getNeighbors() {
    const size_t numPoints = points.size();
    std::vector<std::vector<size_t>> neighbors(numPoints);

    for (size_t i = 0; i < numPoints; ++i) {
      const IdPoint &point = points[i];
      size_t cellIndex = calculateCellIndex(point);
      std::set<size_t> uniqueNeighborIndices;

      const std::vector<size_t> &neighborIndices = neighborCellIndices[cellIndex];
      for (size_t neighborCellIndex : neighborIndices) {
        for (size_t neighborIndex : grid[neighborCellIndex]) {
          if (neighborIndex <= i)
            continue;
          uniqueNeighborIndices.insert(neighborIndex);
        }
      }

      // neighbors[i].reserve(uniqueNeighborIndices.size());
      for (size_t neighborIndex : uniqueNeighborIndices) {
        neighbors[i].push_back(neighborIndex);
      }
    }

    return neighbors;
  }

private:
  /// Calculates the cell index for a given point.
  ///
  /// @param point The point.
  /// @return The cell index.
  ///
  size_t calculateCellIndex(const IdPoint &point) const {
    size_t ix = (size_t)floor(point.x * gridSizeDouble);
    size_t iy = (size_t)floor(point.y * gridSizeDouble);
    size_t iz = (size_t)floor(point.z * gridSizeDouble);

    return iz * gridSizeSquared + iy * gridSize + ix;
  }

  std::vector<size_t> calculateNeighborCellIndices(size_t cellIndex) const {
    const size_t z = cellIndex / gridSizeSquared;
    const size_t y = (cellIndex % gridSizeSquared) / gridSize;
    const size_t x = cellIndex % gridSize;

    std::vector<size_t> neighborIndices;
    // neighborIndices.reserve(27);

    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {

          int nx = (int)x + dx;
          if (nx < 0) {
            nx = gridSizeInt - 1;
          } else if (nx >= gridSizeInt) {
            nx = 0;
          }

          int ny = (int)y + dy;
          if (ny < 0) {
            ny = gridSizeInt - 1;
          } else if (ny >= gridSizeInt) {
            ny = 0;
          }

          int nz = (int)z + dz;
          if (nz < 0) {
            nz = gridSizeInt - 1;
          } else if (nz >= gridSizeInt) {
            nz = 0;
          }

          size_t neighborIndex =
              static_cast<size_t>(nz) * gridSizeSquared + static_cast<size_t>(ny) * gridSize + static_cast<size_t>(nx);
          neighborIndices.push_back(neighborIndex);
        }
      }
    }

    return neighborIndices;
  }

  /// Calculates the distance between two points, accounting for periodicity.
  ///
  /// @param point1 The first point.
  /// @param point2 The second point.
  /// @return The distance between the points.
  ///
  double calculateSqrDistance(const IdPoint &point1, const IdPoint &point2) const {
    double sx = point2.x - point1.x;
    double sy = point2.y - point1.y;
    double sz = point2.z - point1.z;

    sx -= floor(sx + 0.5);
    sy -= floor(sy + 0.5);
    sz -= floor(sz + 0.5);

    return sx * sx + sy * sy + sz * sz;
  }

  ///
  /// Builds the grid and assigns points to grid cells.
  ///
  void buildGrid() {
    grid.resize(numCells);

    neighborCellIndices.resize(numCells);
    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
      neighborCellIndices[cellIndex] = calculateNeighborCellIndices(cellIndex);
    }

    for (size_t i = 0; i < points.size(); ++i) {
      const IdPoint &point = points[i];
      size_t cellIndex = calculateCellIndex(point);
      grid[cellIndex].push_back(i);
    }
  }

  std::vector<IdPoint> points;
  size_t gridSize;
  size_t numCells;
  double gridSizeDouble;
  int gridSizeInt;
  size_t gridSizeSquared;
  std::vector<std::vector<size_t>> grid;
  std::vector<std::vector<size_t>> neighborCellIndices;
};

#endif /* end of include guard: PERIODIC_NN_HPP */

#if 0

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

int main() {

  std::vector<IdPoint> points;
  const size_t numPoints = 4 * 4 * 4;
  size_t i = 0;
  for (size_t iz = 0; iz < 4; ++iz) {
    for (size_t iy = 0; iy < 4; ++iy) {
      for (size_t ix = 0; ix < 4; ++ix) {
        double x = 0.2 + 0.2 * ix;
        double y = 0.2 + 0.2 * iy;
        double z = 0.2 + 0.2 * iz;
        points.push_back({i, x, y, z});
        i++;
      }
    }
  }

  /*
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dist(0.0, 1.0);

std::vector<IdPoint> points;
const size_t numPoints = 20;
for (size_t i = 0; i < numPoints; ++i) {
double x = dist(gen);
double y = dist(gen);
double z = dist(gen);
points.push_back({i, x, y, z});
}
  */

  /*
  std::vector<IdPoint> points;
  const size_t numPoints = 4;
  points.push_back({0, 0.1, 0.1, 0.1});
  points.push_back({1, 0.9, 0.9, 0.1});
  points.push_back({2, 0.9, 0.1, 0.1});
  points.push_back({3, 0.1, 0.9, 0.1});
  */

  // Determine the gridSize based on the density of points
  double averageDistance = std::pow(1.0 / numPoints, 1.0 / 3.0); // should be in ]0 0.5[
  std::cout << "averageDistance: " << averageDistance << std::endl;

  size_t gridSize = std::max<size_t>(static_cast<size_t>(std::ceil(1.0 / averageDistance)), 3);
  // size_t gridSize = 3;
  std::cout << "gridSize: " << gridSize << std::endl;

  PeriodicNearestNeighbors perioNN(points, gridSize);

  auto start = std::chrono::high_resolution_clock::now(); 
  
	std::vector<std::vector<size_t>> neighbors = perioNN.getNeighbors(1.7 * averageDistance);

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Execution time: " << duration << " microseconds" << std::endl;

  for (size_t i = 0; i < neighbors.size(); ++i) {
    std::cout << "Neighbors of point " << i << ": ";
    for (size_t neighbor : neighbors[i]) {
      std::cout << neighbor << " ";
    }
    std::cout << std::endl;
  }

  std::ofstream file("data2.txt");
  for (size_t i = 0; i < neighbors.size(); ++i) {
    for (size_t j : neighbors[i]) {
      file << points[i].x << " " << points[i].y << " " << points[i].z << "\n";
      file << points[j].x << " " << points[j].y << " " << points[j].z << "\n\n";
    }
  }

  return 0;
}

#endif

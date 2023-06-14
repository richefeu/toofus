#ifndef PERIODIC_NN_HPP
#define PERIODIC_NN_HPP

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

struct IdPoint {
  size_t id;
  double x, y, z;
};

class PeriodicNearestNeighbors {
public:
  /**
   * Constructs a PeriodicNearestNeighbors object with the given points and gridSize.
   *
   * @param points The vector of points.
   * @param gridSize The size of the grid.
   */
  PeriodicNearestNeighbors(const std::vector<IdPoint> &t_points, size_t t_gridSize)
      : points(t_points), gridSize(t_gridSize), numCells(t_gridSize * t_gridSize * t_gridSize) {
    buildGrid();
  }

  /**
   * Returns the neighbors of each point within the given radius.
   *
   * @param dmax The search radius.
   * @return A vector of vectors containing the indices of neighbors for each point.
   */
  std::vector<std::vector<size_t>> getNeighbors(double dmax) {
    const size_t numPoints = points.size();
    std::vector<std::vector<size_t>> neighbors(numPoints);

    for (size_t i = 0; i < numPoints; ++i) {
      const IdPoint &point = points[i];
      size_t cellIndex = calculateCellIndex(point);

      std::set<size_t> uniqueNeighborIndices;

      for (size_t neighborCellIndex : getNeighborCellIndices(cellIndex)) {
        for (size_t neighborIndex : grid[neighborCellIndex]) {
          if (neighborIndex <= i) // Skip points with smaller or equal identifiers
            continue;

          const IdPoint &neighbor = points[neighborIndex];
          double sqrDistance = calculateSqrDistance(point, neighbor);

          if (sqrDistance <= dmax * dmax) {
            uniqueNeighborIndices.insert(neighborIndex);
          }
        }
      }

      neighbors[i].reserve(uniqueNeighborIndices.size());
      for (size_t neighborIndex : uniqueNeighborIndices) {
        neighbors[i].push_back(neighborIndex);
      }
    }

    return neighbors;
  }

private:
  /**
   * Calculates the cell index for a given point.
   *
   * @param point The point.
   * @return The cell index.
   */
  size_t calculateCellIndex(const IdPoint &point) const {
double gridSizeDouble = static_cast<double>(gridSize);
    size_t ix = (size_t)floor(point.x * gridSizeDouble);
    size_t iy = (size_t)floor(point.y * gridSizeDouble);
    size_t iz = (size_t)floor(point.z * gridSizeDouble);

    return iz * gridSize * gridSize + iy * gridSize + ix;
  }

  /**
   * Returns the indices of the neighboring cells for a given cell index.
   *
   * @param cellIndex The cell index.
   * @return A vector of neighboring cell indices.
   */
  std::vector<size_t> getNeighborCellIndices(size_t cellIndex) const {

    size_t gridSizeSquared = gridSize * gridSize;
    size_t z = cellIndex / gridSizeSquared;
    size_t y = (cellIndex % gridSizeSquared) / gridSize;
    size_t x = cellIndex % gridSize;

    std::vector<size_t> neighborIndices;
    neighborIndices.reserve(27);

    int gridSizeInt = static_cast<int>(gridSize);
    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {

          int nx = (int)x + dx;
          if (nx < 0)
            nx = gridSizeInt - 1;
          else if (nx >= gridSizeInt)
            nx = 0;

          int ny = (int)y + dy;
          if (ny < 0)
            ny = gridSizeInt - 1;
          else if (ny >= gridSizeInt)
            ny = 0;

          int nz = (int)z + dz;
          if (nz < 0)
            nz = gridSizeInt - 1;
          else if (nz >= gridSizeInt)
            nz = 0;

          size_t neighborIndex =
              static_cast<size_t>(nz) * gridSizeSquared + static_cast<size_t>(ny) * gridSize + static_cast<size_t>(nx);
          neighborIndices.push_back(neighborIndex);
        }
      }
    }

    return neighborIndices;
  }

  /**
   * Calculates the distance between two points, accounting for periodicity.
   *
   * @param point1 The first point.
   * @param point2 The second point.
   * @return The distance between the points.
   */
  double calculateSqrDistance(const IdPoint &point1, const IdPoint &point2) const {
    double sx = point2.x - point1.x;
    double sy = point2.y - point1.y;
    double sz = point2.z - point1.z;

    sx -= floor(sx + 0.5);
    sy -= floor(sy + 0.5);
    sz -= floor(sz + 0.5);

    return sx * sx + sy * sy + sz * sz;
  }

  /**
   * Builds the grid and assigns points to grid cells.
   */
  void buildGrid() {
    grid.resize(numCells);

    for (size_t i = 0; i < points.size(); ++i) {
      const IdPoint &point = points[i];
      size_t cellIndex = calculateCellIndex(point);
      grid[cellIndex].push_back(i);
    }
  }

  std::vector<IdPoint> points;
  size_t gridSize;
  size_t numCells;
  std::vector<std::vector<size_t>> grid;
};

#endif /* end of include guard: PERIODIC_NN_HPP */

#if 0

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

  std::vector<std::vector<size_t>> neighbors = perioNN.getNeighbors(1.7 * averageDistance);

  for (size_t i = 0; i < neighbors.size(); ++i) {
    std::cout << "Neighbors of point " << i << ": ";
    for (size_t neighbor : neighbors[i]) {
      std::cout << neighbor << " ";
    }
    std::cout << std::endl;
  }

	std::ofstream file ("data.txt");
  for (size_t i = 0; i < neighbors.size(); ++i) {
    for (size_t j : neighbors[i]) {
      file << points[i].x << " " << points[i].y << " " << points[i].z << "\n";
			file << points[j].x << " " << points[j].y << " " << points[j].z << "\n\n";
    }
  }



  return 0;
}

#endif

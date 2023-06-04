#ifndef OBB_2D_QUADTREE_HPP
#define OBB_2D_QUADTREE_HPP

#include <array>
#include <memory>
#include <vector>

#include "OBB_2D.hpp"

class OBB_2D_quadtree {
private:
  static constexpr int MAX_CAPACITY = 4;

  struct TreeNode {
    OBB_2D bounds;
    std::vector<OBB_2D> objects;
    std::array<std::unique_ptr<TreeNode>, 4> children;

    TreeNode(const OBB_2D &bbox) : bounds(bbox), objects() {}
  };

  std::unique_ptr<TreeNode> root;

  void insert(const OBB_2D &object, std::unique_ptr<TreeNode> &node) {
    if (!node->bounds.intersect(object))
      return;

    if (node->objects.size() < MAX_CAPACITY) {
      node->objects.push_back(object);
      return;
    }

    if (node->children[0] == nullptr) {
      splitNode(node);
    }

    for (auto &child : node->children) {
      insert(object, child);
    }
  }

  void splitNode(std::unique_ptr<TreeNode> &node) {

    const double subWidth = node->bounds.width / 2.0;
    const double subHeight = node->bounds.height / 2.0;
    const double xMid = node->bounds.position.x;
    const double yMid = node->bounds.position.y;

    std::array<OBB_2D, 4> childBounds = {
        OBB_2D(xMid - subWidth / 2.0, yMid - subHeight / 2.0, subWidth, subHeight, 0.0),
        OBB_2D(xMid + subWidth / 2.0, yMid - subHeight / 2.0, subWidth, subHeight, 0.0),
        OBB_2D(xMid - subWidth / 2.0, yMid + subHeight / 2.0, subWidth, subHeight, 0.0),
        OBB_2D(xMid + subWidth / 2.0, yMid + subHeight / 2.0, subWidth, subHeight, 0.0)};

    for (size_t i = 0; i < 4; ++i) {
      node->children[i] = std::make_unique<TreeNode>(std::move(childBounds[i]));
    }

    for (const auto &obj : node->objects) {
      for (auto &child : node->children) {
        insert(obj, child);
      }
    }

    node->objects.clear();
  }

  void query(const OBB_2D &range, const std::unique_ptr<TreeNode> &node, std::vector<OBB_2D> &result) const {
    if (!node->bounds.intersect(range))
      return;

    // Reserve capacity for the result vector to avoid frequent reallocations
    result.reserve(result.size() + node->objects.size());

    for (const auto &obj : node->objects) {
      if (obj.intersect(range))
        result.push_back(obj);
    }

    for (const auto &child : node->children) {
      if (child != nullptr) {
        query(range, child, result);
      }
    }
  }

public:
  OBB_2D_quadtree(const OBB_2D &bounds) : root(std::make_unique<TreeNode>(bounds)) {}

  void insert(const OBB_2D &object) { insert(object, root); }

  std::vector<OBB_2D> query(const OBB_2D &range) const {
    std::vector<OBB_2D> result;
    query(range, root, result);
    return result;
  }
};

#endif /* end of include guard: OBB_2D_QUADTREE_HPP */

#if 0

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

void generateSVG(const std::vector<OBB_2D> &obbsBlack, const std::vector<OBB_2D> &obbsRed, double xMin, double yMin,
                 double xMax, double yMax) {
  double scaleFactor = 800.0 / std::max(xMax - xMin, yMax - yMin);

  std::ofstream file("output.svg");
  if (!file.is_open()) {
    std::cout << "Failed to open output.svg file" << std::endl;
    return;
  }

  // Calculate the adjusted coordinates of the OBBs
  auto adjustCoordinates = [&](double x, double y) {
    double adjustedX = (x - xMin) * scaleFactor;
    double adjustedY = ((yMax - yMin) - (y - yMin)) * scaleFactor;
    return std::make_pair(adjustedX, adjustedY);
  };

  // SVG header
  file << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"800\" height=\"800\">" << std::endl;

  // Draw the space limits
  double adjustedXMin, adjustedYMin, adjustedWidth, adjustedHeight;
  std::tie(adjustedXMin, adjustedYMin) = adjustCoordinates(xMin, yMin);
  std::tie(adjustedWidth, adjustedHeight) = adjustCoordinates(xMax, yMax);
  file << "<rect x=\"" << adjustedXMin << "\" y=\"" << adjustedYMin << "\" width=\"" << adjustedWidth << "\" height=\""
       << adjustedHeight << "\" fill=\"none\" stroke=\"black\" />" << std::endl;

  // Draw OBBs in black color
  for (const OBB_2D &obb : obbsBlack) {
    double adjustedX, adjustedY;
    std::tie(adjustedX, adjustedY) = adjustCoordinates(obb.position.x, obb.position.y);
    double adjustedWidth = obb.width * scaleFactor;
    double adjustedHeight = obb.height * scaleFactor;
    double adjustedAngle = obb.rotation;

    file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2
         << "\" width=\"" << adjustedWidth << "\" height=\"" << adjustedHeight << "\" transform=\"rotate("
         << adjustedAngle * 180.0 / M_PI << " " << adjustedX << " " << adjustedY << ")\" fill=\"black\" />"
         << std::endl;
  }

  // Draw OBBs in red color
  for (const OBB_2D &obb : obbsRed) {
    double adjustedX, adjustedY;
    std::tie(adjustedX, adjustedY) = adjustCoordinates(obb.position.x, obb.position.y);
    double adjustedWidth = obb.width * scaleFactor;
    double adjustedHeight = obb.height * scaleFactor;
    double adjustedAngle = obb.rotation;

    file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2
         << "\" width=\"" << adjustedWidth << "\" height=\"" << adjustedHeight << "\" transform=\"rotate("
         << adjustedAngle * 180.0 / M_PI << " " << adjustedX << " " << adjustedY << ")\" fill=\"red\" />" << std::endl;
  }

  // SVG footer
  file << "</svg>" << std::endl;

  file.close();
  std::cout << "SVG file generated: output.svg" << std::endl;
}

int main() {
  // Create a AABB_2D_quadtree with a bounding box
  OBB_2D bounds(0, 0, 16, 16, 0);
  OBB_2D_quadtree quadTree(bounds);

  // Set up random number generation
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> distribution(-8.0, 8.0);
  std::uniform_real_distribution<double> angleDistribution(0.0, M_PI);

  std::vector<OBB_2D> obbs;

  for (double centerX = -8; centerX <= 8; centerX += 16. / 25.) {
    for (double centerY = -8; centerY <= 8; centerY += 16. / 25.) {
      double width = 0.8;
      double height = 0.4;
      double angle = angleDistribution(gen);

      OBB_2D obb(centerX, centerY, width, height, angle);
      quadTree.insert(obb);
      obbs.push_back(obb);
    }
  }

  // Query the QuadTree for intersecting OBBs
  OBB_2D queryRange(2, 0, 12, 1.5, 0.78);
  std::vector<OBB_2D> intersectingOBBs = quadTree.query(queryRange);

  // Generate SVG file with colored OBBs
  generateSVG(obbs, intersectingOBBs, -8, -8, 8, 8);

  return 0;
}

#endif

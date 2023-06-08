#ifndef AABB_2D_QUADTREE_HPP
#define AABB_2D_QUADTREE_HPP

#include <array>
#include <memory>
#include <stack>
#include <vector>

#include "AABB_2D.hpp"

class AABB_2D_quadtree {
public:
  AABB_2D_quadtree(const AABB_2D &bounds, int maxCapacity = 4, int maxDepth = 10)
      : MAX_CAPACITY(maxCapacity), MAX_DEPTH(maxDepth), root(std::make_unique<TreeNode>(bounds)) {}

  void insert(const AABB_2D &object) { insert(object, root, 0); }

  std::vector<AABB_2D> query(const AABB_2D &range) const {
    std::vector<AABB_2D> result;
    query(range, root, result);
    return result;
  }

  void setMaxCapacity(int maxCapacity) { MAX_CAPACITY = maxCapacity; }

  void setMaxDepth(int maxDepth) { MAX_DEPTH = maxDepth; }

private:
  int MAX_CAPACITY;
  int MAX_DEPTH;

  struct TreeNode {
    AABB_2D bounds;
    std::vector<AABB_2D> objects;
    std::array<std::unique_ptr<TreeNode>, 4> children;

    TreeNode(const AABB_2D &bbox) : bounds(bbox), objects() {}
  };

  std::unique_ptr<TreeNode> root;

  void insert(const AABB_2D &object, std::unique_ptr<TreeNode> &node, int depth) {
    if (depth >= MAX_DEPTH || !node->bounds.intersect(object))
      return;

    if (node->objects.size() < MAX_CAPACITY) {
      node->objects.push_back(object);
      return;
    }

    if (node->children[0] == nullptr) {
      splitNode(node);
    }

    for (auto &child : node->children) {
      insert(object, child, depth + 1);
    }
  }

  void splitNode(std::unique_ptr<TreeNode> &node) {
    const double xMid = (node->bounds.min.x + node->bounds.max.x) / 2;
    const double yMid = (node->bounds.min.y + node->bounds.max.x) / 2;

    node->children[0] = std::make_unique<TreeNode>(AABB_2D(node->bounds.min.x, xMid, node->bounds.min.y, yMid));
    node->children[1] = std::make_unique<TreeNode>(AABB_2D(xMid, node->bounds.max.x, node->bounds.min.y, yMid));
    node->children[2] = std::make_unique<TreeNode>(AABB_2D(node->bounds.min.x, xMid, yMid, node->bounds.max.y));
    node->children[3] = std::make_unique<TreeNode>(AABB_2D(xMid, node->bounds.max.x, yMid, node->bounds.max.y));

    for (const auto &obj : node->objects) {
      for (auto &child : node->children) {
        insert(obj, child, 0);
      }
    }

    node->objects.clear();
  }

  void query(const AABB_2D &range, const std::unique_ptr<TreeNode> &node, std::vector<AABB_2D> &result) const {

    if (!node || !node->bounds.intersect(range))
      return;

    for (const auto &obj : node->objects) {
      if (obj.intersect(range))
        result.push_back(obj);
    }

    for (const auto &child : node->children) {
      query(range, child, result);
    }
  }
};

#endif /* end of include guard: AABB_2D_QUADTREE_HPP */

#if 0

#include "Mth.hpp"
#include <fstream>
#include <iostream>

void generateSVG(const std::vector<AABB_2D> &aabbsBlack, const std::vector<AABB_2D> &aabbsRed, double xMin, double yMin,
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
  double adjustedXMin, adjustedYMin, adjustedXMax, adjustedYMax;
  std::tie(adjustedXMin, adjustedYMin) = adjustCoordinates(xMin, yMin);
  std::tie(adjustedXMax, adjustedYMax) = adjustCoordinates(xMax, yMax);
  file << "<rect x=\"" << adjustedXMin << "\" y=\"" << adjustedYMin << "\" width=\"" << adjustedXMax - adjustedXMin
       << "\" height=\"" << adjustedYMax - adjustedYMin << "\" fill=\"none\" stroke=\"black\" />" << std::endl;

  // Draw OBBs in black color
  for (const AABB_2D &aabb : aabbsBlack) {
    double adjustedX, adjustedY;
    std::tie(adjustedX, adjustedY) =
        adjustCoordinates(0.5 * (aabb.min.x + aabb.max.x), 0.5 * (aabb.min.y + aabb.max.y));
    double adjustedWidth = (aabb.max.x - aabb.min.x) * scaleFactor;
    double adjustedHeight = (aabb.max.y - aabb.min.y) * scaleFactor;
    // double adjustedAngle = obb.rotation;

    file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2
         << "\" width=\"" << adjustedWidth << "\" height=\"" << adjustedHeight << "\" fill=\"none\" stroke=\"black\"/>"
         << std::endl;
  }

  // Draw OBBs in red color
  for (const AABB_2D &aabb : aabbsRed) {
    double adjustedX, adjustedY;
    std::tie(adjustedX, adjustedY) =
        adjustCoordinates(0.5 * (aabb.min.x + aabb.max.x), 0.5 * (aabb.min.y + aabb.max.y));
    double adjustedWidth = (aabb.max.x - aabb.min.x) * scaleFactor;
    double adjustedHeight = (aabb.max.y - aabb.min.y) * scaleFactor;
    // double adjustedAngle = obb.rotation;

    file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2
         << "\" width=\"" << adjustedWidth << "\" height=\"" << adjustedHeight << "\" fill=\"none\" stroke=\"red\" />"
         << std::endl;
  }

  // SVG footer
  file << "</svg>" << std::endl;

  file.close();
  std::cout << "SVG file generated: output.svg" << std::endl;
}

int main() {
  // Create a AABB_2D_quadtree with a bounding box
  AABB_2D bounds(-10, 10, -10, 10);
  AABB_2D_quadtree quadTree(bounds, 16, 10);

  std::vector<AABB_2D> aabbs;

  for (double centerX = -8; centerX <= 8; centerX += .5) {
    for (double centerY = -8; centerY <= 8; centerY += .5) {
      double width = Mth::random(0.5, .6);
      double height = Mth::random(0.5, 0.6);
      double dx = Mth::random(-0.2, 0.2);
      double dy = Mth::random(-0.2, 0.2);

      AABB_2D aabb(dx + centerX - width / 2.0, dx + centerX + width / 2.0, dy + centerY - height / 2.0,
                   dy + centerY + height / 2.0);
      quadTree.insert(aabb);
      aabbs.push_back(aabb);
    }
  }

  // Query the AABB_2D_quadtree for intersecting AABBs
  AABB_2D queryRange(-5, 1, -4.5, -2);
  std::vector<AABB_2D> intersectingAABBs = quadTree.query(queryRange);

  generateSVG(aabbs, intersectingAABBs, -8, -8, 8, 8);

  return 0;
}

#endif

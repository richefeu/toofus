#ifndef AABB_2D_QUADTREE_HPP
#define AABB_2D_QUADTREE_HPP

#include <memory>
#include <vector>
#include <array>

#include "AABB_2D.hpp"

class AABB_2D_quadtree {
private:
  static constexpr int MAX_CAPACITY = 4;

  struct TreeNode {
    AABB_2D bounds;
    std::vector<AABB_2D> objects;
    std::array<std::unique_ptr<TreeNode>, 4> children;

    TreeNode(const AABB_2D &bbox) : bounds(bbox), objects() {}
  };

  std::unique_ptr<TreeNode> root;

  void insert(const AABB_2D &object, std::unique_ptr<TreeNode> &node) {
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
    const double xMid = (node->bounds.min.x + node->bounds.max.x) / 2;
    const double yMid = (node->bounds.min.y + node->bounds.max.x) / 2;

    node->children[0] = std::make_unique<TreeNode>(AABB_2D(node->bounds.min.x, xMid, node->bounds.min.y, yMid));
    node->children[1] = std::make_unique<TreeNode>(AABB_2D(xMid, node->bounds.max.x, node->bounds.min.y, yMid));
    node->children[2] = std::make_unique<TreeNode>(AABB_2D(node->bounds.min.x, xMid, yMid, node->bounds.max.y));
    node->children[3] = std::make_unique<TreeNode>(AABB_2D(xMid, node->bounds.max.x, yMid, node->bounds.max.y));

    for (const auto &obj : node->objects) {
      for (auto &child : node->children) {
        insert(obj, child);
      }
    }

    node->objects.clear();
  }

  void query(const AABB_2D &range, const std::unique_ptr<TreeNode> &node, std::vector<AABB_2D> &result) const {
    if (!node->bounds.intersect(range))
      return;

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
  AABB_2D_quadtree(const AABB_2D &bounds) : root(std::make_unique<TreeNode>(bounds)) {}

  void insert(const AABB_2D &object) { insert(object, root); }

  std::vector<AABB_2D> query(const AABB_2D &range) const {
    std::vector<AABB_2D> result;
    query(range, root, result);
    return result;
  }
};

#endif /* end of include guard: AABB_2D_QUADTREE_HPP */

#if 0

#include <iostream>

int main() {
  // Create a AABB_2D_quadtree with a bounding box
  AABB_2D bounds(-10, 10, -10, 10);
  AABB_2D_quadtree quadTree(bounds);

  // Insert some AABBs into the AABB_2D_quadtree
  quadTree.insert(AABB_2D(-5, -3, -3, -1));
  quadTree.insert(AABB_2D(2, 5, 2, 4));
  quadTree.insert(AABB_2D(0, 1, 0, 1));
  quadTree.insert(AABB_2D(-8, -6, -6, -4));

  // Query the AABB_2D_quadtree for intersecting AABBs
  AABB_2D queryRange(-7, -4, -4.5, 0);
  std::vector<AABB_2D> intersectingAABBs = quadTree.query(queryRange);

  // Print the intersecting AABBs
  for (const AABB_2D &aabb : intersectingAABBs) {
    std::cout << "Intersecting AABB: (" << aabb.min.x << ", " << aabb.max.x << ", " << aabb.min.y << ", " << aabb.max.y
              << ")" << std::endl;
  }

  return 0;
}

#endif

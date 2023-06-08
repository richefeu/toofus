#ifndef OBB_2D_QUADTREE_HPP
#define OBB_2D_QUADTREE_HPP

#include <array>
#include <memory>
#include <vector>

#include "OBB_2D.hpp"

class OBB_2D_quadtree {
public:
  OBB_2D_quadtree(const OBB_2D &bounds, int maxCapacity = 4, int maxDepth = 10)
      : MAX_CAPACITY(maxCapacity), MAX_DEPTH(maxDepth), root(std::make_unique<TreeNode>(bounds)) {}

  void insert(const OBB_2D &object, size_t objectId) { insert(object, objectId, root, 0); }

  std::vector<OBB_2D> query(const OBB_2D &range) const {
    std::vector<OBB_2D> result;
    query(range, root, result);
    return result;
  }

  std::vector<OBB_2D> query(size_t objectId) const {
    std::vector<OBB_2D> result;

    if (objectId < objectInfo.size()) {
      OBB_2D range = objectInfo[objectId].node->objects[objectInfo[objectId].locObjectId];
      query(range, root, result, objectId);
    }

    return result;
  }

  std::vector<size_t> queryIDs(size_t objectId) const {
    std::vector<size_t> result;

    if (objectId < objectInfo.size()) {
      OBB_2D range = objectInfo[objectId].node->objects[objectInfo[objectId].locObjectId];
      queryIDs(range, root, result, objectId);
    }

    return result;
  }
	
  std::vector<size_t> queryIDsFrom(size_t objectId) const {
    std::vector<size_t> result;

    if (objectId < objectInfo.size()) {
      OBB_2D range = objectInfo[objectId].node->objects[objectInfo[objectId].locObjectId];
      queryIDsFrom(range, root, result, objectId);
    }

    return result;
  }

  // Setting maximum capacity and maximum depth is a matter of optimisation
  void setMaxCapacity(int maxCapacity) { MAX_CAPACITY = maxCapacity; }
  void setMaxDepth(int maxDepth) { MAX_DEPTH = maxDepth; }

  // Set the size of the objectInfo vector
  // It is better to set this size before feeling the quadTree when this size is known
  void setObjectInfoSize(size_t size) { objectInfo.resize(size); }

private:
  int MAX_CAPACITY;
  int MAX_DEPTH;

  struct TreeNode {
    OBB_2D bounds;
    std::vector<OBB_2D> objects;
    std::vector<size_t> objectIds; // Stores the object IDs
    std::array<std::unique_ptr<TreeNode>, 4> children;

    TreeNode(const OBB_2D &bbox) : bounds(bbox), objects() {}
  };

  std::unique_ptr<TreeNode> root;

  struct objPos {
    TreeNode *node;
    size_t locObjectId;

    objPos(TreeNode *t_node, size_t t_locId) : node(t_node), locObjectId(t_locId) {}
    objPos() : node(nullptr), locObjectId(0) {}
  };

  std::vector<objPos> objectInfo;

  void insert(const OBB_2D &object, size_t objectId, std::unique_ptr<TreeNode> &node, int depth) {
    if (depth >= MAX_DEPTH || !node->bounds.contains(object))
      return;

    if (node->objects.size() < MAX_CAPACITY) {
      node->objects.push_back(object);
      node->objectIds.push_back(objectId);
      if (objectId >= objectInfo.size()) {
        objectInfo.resize(objectId + 1);
      }
      objectInfo[objectId] = objPos(node.get(), node->objects.size() - 1);
      return;
    }

    if (node->children[0] == nullptr) {
      splitNode(node);
    }

    for (auto &child : node->children) {
      insert(object, objectId, child, depth + 1);
    }
  }

  void splitNode(std::unique_ptr<TreeNode> &node) {

    const double subWidth = node->bounds.width / 2.0;
    const double subHeight = node->bounds.height / 2.0;
    const double halfSubWidth = subWidth / 2.0;
    const double halfSubHeight = subHeight / 2.0;
    const double xMid = node->bounds.position.x;
    const double yMid = node->bounds.position.y;
    const double parentRotation = node->bounds.rotation;

    const double c = cos(parentRotation);
    const double s = sin(parentRotation);

    const double Wx = c * halfSubWidth;
    const double Wy = s * halfSubWidth;
    const double Hx = -s * halfSubWidth;
    const double Hy = c * halfSubWidth;

    std::array<OBB_2D, 4> childBounds = {OBB_2D(xMid - Wx - Hx, yMid - Wy - Hy, subWidth, subHeight, parentRotation),
                                         OBB_2D(xMid + Wx - Hx, yMid + Wy - Hy, subWidth, subHeight, parentRotation),
                                         OBB_2D(xMid - Wx + Hx, yMid - Wy + Hy, subWidth, subHeight, parentRotation),
                                         OBB_2D(xMid + Wx + Hx, yMid + Wy + Hy, subWidth, subHeight, parentRotation)};

    node->children[0] = std::make_unique<TreeNode>(std::move(childBounds[0]));
    node->children[1] = std::make_unique<TreeNode>(std::move(childBounds[1]));
    node->children[2] = std::make_unique<TreeNode>(std::move(childBounds[2]));
    node->children[3] = std::make_unique<TreeNode>(std::move(childBounds[3]));

    for (size_t i = 0; i < node->objects.size(); ++i) {
      const auto &obj = node->objects[i];
      const size_t objectId = node->objectIds[i];
      for (auto &child : node->children) {
        insert(obj, objectId, child, 0);
      }
    }

    node->objects.clear();
    node->objectIds.clear();
  }

  void query(const OBB_2D &range, const std::unique_ptr<TreeNode> &node, std::vector<OBB_2D> &result,
             int fromId = -1) const {
    if (!node || !node->bounds.intersects(range))
      return;

    if (fromId < 0) {

      for (size_t i = 0; i < node->objects.size(); ++i) {
        const auto &obj = node->objects[i];
        if (obj.intersects(range)) {
          result.push_back(obj);
        }
      }

    } else {

      size_t Id = static_cast<size_t>(std::abs(fromId));

      for (size_t i = 0; i < node->objects.size(); ++i) {
        const auto &obj = node->objects[i];
        const size_t objectId = node->objectIds[i];

        if (objectId > Id && obj.intersects(range)) {
          result.push_back(obj);
        }
      }
    }

    for (const auto &child : node->children) {
      query(range, child, result, fromId);
    }
  }

  void queryIDs(const OBB_2D &range, const std::unique_ptr<TreeNode> &node, std::vector<size_t> &result,
                size_t theId) const {
    if (!node || !node->bounds.intersects(range))
      return;

    for (size_t i = 0; i < node->objects.size(); ++i) {
      const auto &obj = node->objects[i];
      const size_t objectId = node->objectIds[i];

      if (objectId != theId && obj.intersects(range)) {
        result.push_back(objectId);
      }
    }

    for (const auto &child : node->children) {
      queryIDs(range, child, result, theId);
    }
  }
	
  void queryIDsFrom(const OBB_2D &range, const std::unique_ptr<TreeNode> &node, std::vector<size_t> &result,
                size_t fromId) const {
    if (!node || !node->bounds.intersects(range))
      return;

    for (size_t i = 0; i < node->objects.size(); ++i) {
      const auto &obj = node->objects[i];
      const size_t objectId = node->objectIds[i];

      if (objectId > fromId && obj.intersects(range)) {
        result.push_back(objectId);
      }
    }

    for (const auto &child : node->children) {
      queryIDsFrom(range, child, result, fromId);
    }
  }
};

#endif /* end of include guard: OBB_2D_QUADTREE_HPP */

#if 1

#include <chrono>
#include <fstream>
#include <iostream>

#include "Mth.hpp"
#include "message.hpp"

void generateSVG(const char *fname, OBB_2D &obbBlue, const std::vector<OBB_2D> &obbsBlack,
                 const std::vector<OBB_2D> &obbsRed, double xMin, double yMin, double xMax, double yMax) {
  double scaleFactor = 800.0 / std::max(xMax - xMin, yMax - yMin);

  std::ofstream file(fname);
  if (!file.is_open()) {
    std::cout << "Failed to open " << fname << " file" << std::endl;
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

  // Draw the blue OBB
  double adjustedX, adjustedY;
  std::tie(adjustedX, adjustedY) = adjustCoordinates(obbBlue.position.x, obbBlue.position.y);
  double adjustedWidth = obbBlue.width * scaleFactor;
  double adjustedHeight = obbBlue.height * scaleFactor;
  double adjustedAngle = obbBlue.rotation;
  file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2 << "\" width=\""
       << adjustedWidth << "\" height=\"" << adjustedHeight << "\" transform=\"rotate(" << -adjustedAngle * 180.0 / M_PI
       << " " << adjustedX << " " << adjustedY << ")\" fill=\"blue\" stroke=\"blue\" />" << std::endl;

  // Draw OBBs in black color
  for (const OBB_2D &obb : obbsBlack) {
    double adjustedX, adjustedY;
    std::tie(adjustedX, adjustedY) = adjustCoordinates(obb.position.x, obb.position.y);
    double adjustedWidth = obb.width * scaleFactor;
    double adjustedHeight = obb.height * scaleFactor;
    double adjustedAngle = obb.rotation;

    file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2
         << "\" width=\"" << adjustedWidth << "\" height=\"" << adjustedHeight << "\" transform=\"rotate("
         << -adjustedAngle * 180.0 / M_PI << " " << adjustedX << " " << adjustedY
         << ")\" fill=\"none\" stroke=\"black\" />" << std::endl;
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
         << -adjustedAngle * 180.0 / M_PI << " " << adjustedX << " " << adjustedY
         << ")\" fill=\"none\" stroke=\"red\" />" << std::endl;
  }

  // SVG footer
  file << "</svg>" << std::endl;

  file.close();
  std::cout << "SVG file generated: " << fname << std::endl;
}

#include <algorithm>

int main() {

  // Set up random number generation
  std::vector<OBB_2D> obbs;
  for (double centerX = -16; centerX <= 16; centerX += 1.0) {
    for (double centerY = -16; centerY <= 16; centerY += 1.0) {
      double width = Mth::random(0.8, 1.6);
      double height = Mth::random(0.8, 1.6);
      double angle = Mth::random(0.0, M_PI);
      double dx = Mth::random(-0.2, 0.2);
      double dy = Mth::random(-0.2, 0.2);

      OBB_2D obb(dx + centerX, dy + centerY, width, height, angle);

      obbs.push_back(obb);
    }
  }
  std::cout << "nb OBBs: " << obbs.size() << std::endl;
  size_t objectId = 3;//obbs.size() / 2 + 130;

  // O(N^2) complexity
  std::vector<std::pair<size_t, size_t>> VerletListSlow;
  std::vector<OBB_2D> intersOBBs1;
  auto start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < obbs.size(); i++) {
    for (size_t j = i + 1; j < obbs.size(); j++) {
      if (obbs[i].intersects(obbs[j])) {
        VerletListSlow.push_back(std::make_pair(i, j));
        if (i == objectId)
          intersOBBs1.push_back(obbs[j]);
      }
    }
  }
  auto end = std::chrono::steady_clock::now();
  std::cout << "O(N^2) duration: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << " micro-seconds" << std::endl;

  __SHOW(VerletListSlow.size());
  __SHOW(intersOBBs1.size());

  // Generate SVG file with colored OBBs
  generateSVG("withBruteForce.svg", obbs[objectId], obbs, intersOBBs1, -18, -18, 18, 18);

  // Building the quadtree
  start = std::chrono::steady_clock::now();
  OBB_2D bounds(0, 0, 33, 33, 0);
  OBB_2D_quadtree quadTree(bounds, 128, 30);
  quadTree.setObjectInfoSize(obbs.size());
  for (size_t i = 0; i < obbs.size(); i++) {
    quadTree.insert(obbs[i], i);
  }

  std::vector<std::pair<size_t, size_t>> VerletListFast;
  std::vector<OBB_2D> intersectingOBBs;
  for (size_t i = 0; i < obbs.size(); i++) {
    std::vector<size_t> tmp = quadTree.queryIDs(i);
		std::sort(tmp.begin(), tmp.end());
    for (size_t k = 0; k < tmp.size(); k++) {
      VerletListFast.push_back(std::make_pair(i, tmp[k]));
    }

    if (i == objectId) {
      intersectingOBBs = quadTree.query(objectId);
      __SHOW(tmp.size());
			__SHOW(obbs[tmp[tmp.size()-1]].position);
			__SHOW(intersectingOBBs[tmp.size()-1].position);
    }
  }

  end = std::chrono::steady_clock::now();
  std::cout << "quadtree duration: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << " micro-seconds" << std::endl;

  __SHOW(VerletListFast.size());
  __SHOW(intersectingOBBs.size());
	
	for (size_t i = 0 ; i < 20; i++) {
		__SHOW(VerletListSlow[i].first);
		__SHOW(VerletListSlow[i].second);
		__SHOW(VerletListFast[i].first);
		__SHOW(VerletListFast[i].second);
		std::cerr << std::endl;
	}
	

  // Generate SVG file with colored OBBs
  generateSVG("withQuadTree.svg", obbs[objectId], obbs, intersectingOBBs, -18, -18, 18, 18);

  return 0; // ##############################

  // Vary the maximum capacity
  size_t MAX_DEPTH = 20;
  for (size_t MAX_CAPACITY = 1; MAX_CAPACITY <= 50; MAX_CAPACITY++) {
    // Create a new OBB_2D quadtree
    OBB_2D bounds(0, 0, 33, 33, 0);
    OBB_2D_quadtree quadTree(bounds, MAX_CAPACITY, MAX_DEPTH);
    quadTree.setObjectInfoSize(obbs.size());

    // Insert the OBBs into the quadtree
    for (size_t i = 0; i < obbs.size(); i++) {
      quadTree.insert(obbs[i], i);
    }

    // Measure the execution time of a query
    auto start = std::chrono::steady_clock::now();
    std::vector<OBB_2D> queriedObbs = quadTree.query(objectId);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // Print the result
    std::cout << MAX_CAPACITY << ' ' << queriedObbs.size() << ' ' << duration << std::endl;

    /*
std::cout << "Query result for MAX_CAPACITY = " << MAX_CAPACITY << ":" << std::endl;
std::cout << "Nb OBBs found: " << queriedObbs.size() << std::endl;
std::cout << "Execution time: " << duration << " microseconds" << std::endl;
std::cout << std::endl;
    */
  }

  return 0;
}

#endif

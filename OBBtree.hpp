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

#ifndef OBB_TREE_HPP
#define OBB_TREE_HPP

#include <cmath>
//#include <set>
#include <algorithm> // for std::sort
#include <utility>   // for std::pair
#include <vector>

#include "OBB.hpp"
#include "quat.hpp"
#include "vec3.hpp"

struct OBBbundle {
  OBB obb;
  int Id;
  std::vector<vec3r> points; // used for fitting sets of 'leafs'

  OBBbundle(const OBB &obb_, const int id) : obb(obb_), Id(id) {}
};

struct OBBnode {
  OBB boundary;
  OBBnode *first;
  OBBnode *second;
  int refId;
  OBBnode() : boundary(), first(nullptr), second(nullptr), refId(-1) {}
  bool isLeaf() { return (refId >= 0); }
};

class OBBtree {
public:
  OBBnode *root; // the root of the OBB-tree

  OBBtree() : root(nullptr) {}
  ~OBBtree() { reset(root); }

  void reset(OBBnode *node) {
    if (node == nullptr)
      return;

    // first delete subtrees
    reset(node->first);
    reset(node->second);

    // then delete the node
    delete node;
    node = nullptr;
  }

  static void subdivide(std::vector<OBBbundle> &OBBs, std::vector<OBBbundle> &first_OBBs,
                        std::vector<OBBbundle> &second_OBBs) {
    if (OBBs.size() <= 1)
      return;
    if (OBBs.size() == 2) {
      first_OBBs.push_back(OBBs[0]);
      second_OBBs.push_back(OBBs[1]);
      return;
    }

    // find the axis of greatest extension
    mat9r C = OBBtree::getCovarianceMatrix(OBBs);
    mat9r eigvec;
    vec3r eigval;
    C.sorted_sym_eigen(eigvec, eigval);
    vec3r u(eigvec.xx, eigvec.yx, eigvec.zx);

    // project onto this axis
    std::vector<double> proj;
    for (size_t i = 0; i < OBBs.size(); i++) {
      proj.push_back(OBBs[i].obb.center * u);
    }

    // compute median value of projections
    std::vector<double> ordered_proj(proj);
    std::sort(ordered_proj.begin(), ordered_proj.end());
    double median;
    size_t mid = OBBs.size() / 2;
    if (OBBs.size() % 2 == 0) { // even (pair)
      median = 0.5 * (ordered_proj[mid - 1] + ordered_proj[mid]);
    } else { // odd (impair)
      median = ordered_proj[mid];
    }

    // distribute
    for (size_t i = 0; i < proj.size(); i++) {
      if (proj[i] < median) {
        first_OBBs.push_back(OBBs[i]);
      } else {
        second_OBBs.push_back(OBBs[i]);
      }
    }
  }

  // Build the covariance matrix with the points in OBB bundles
  static mat9r getCovarianceMatrix(std::vector<OBBbundle> &OBBs) {
    vec3r mu;
    mat9r C;

    // loop over the points to find the mean point
    // location
    size_t nbP = 0;
    for (size_t io = 0; io < OBBs.size(); io++) {
      for (size_t p = 0; p < OBBs[io].points.size(); p++) {
        mu += OBBs[io].points[p];
        nbP++;
      }
    }

    if (nbP == 0) {
      std::cerr << "@getCovarianceMatrix, no points in the OBBs!!\n";
      return C;
    }
    mu /= (double)nbP;

    // loop over the points again to build the
    // covariance matrix.  Note that we only have
    // to build terms for the upper trianglular
    // portion since the matrix is symmetric
    double cxx = 0.0, cxy = 0.0, cxz = 0.0, cyy = 0.0, cyz = 0.0, czz = 0.0;

    for (size_t io = 0; io < OBBs.size(); io++) {
      for (size_t p = 0; p < OBBs[io].points.size(); p++) {

        vec3r pt = OBBs[io].points[p];
        cxx += pt.x * pt.x - mu.x * mu.x;
        cxy += pt.x * pt.y - mu.x * mu.y;
        cxz += pt.x * pt.z - mu.x * mu.z;
        cyy += pt.y * pt.y - mu.y * mu.y;
        cyz += pt.y * pt.z - mu.y * mu.z;
        czz += pt.z * pt.z - mu.z * mu.z;
      }
    }

    // now build the covariance matrix
    C.xx = cxx;
    C.xy = cxy;
    C.xz = cxz;
    C.yx = cxy;
    C.yy = cyy;
    C.yz = cyz;
    C.zx = cxz;
    C.zy = cyz;
    C.zz = czz;

    return C;
  }

  static OBB fitOBB(std::vector<OBBbundle> &OBBs, double radius = 0.0) {
    OBB fittedObb;
    if (OBBs.empty()) {
      std::cerr << "@fitOBB, OBBs is empty!!\n";
      return fittedObb;
    }

    // compute the covariance matrix
    mat9r C = OBBtree::getCovarianceMatrix(OBBs);

    // ==== set the OBB parameters from the covariance matrix
    // extract the eigenvalues and eigenvectors from C
    mat9r eigvec;
    vec3r eigval;
    C.sym_eigen(eigvec, eigval);

    // find the right, up and forward vectors from the eigenvectors
    vec3r r(eigvec.xx, eigvec.yx, eigvec.zx);
    vec3r u(eigvec.xy, eigvec.yy, eigvec.zy);
    vec3r f(eigvec.xz, eigvec.yz, eigvec.zz);
    r.normalize();
    u.normalize(), f.normalize();

    // now build the bounding box extents in the rotated frame
    vec3r minim(1e20, 1e20, 1e20), maxim(-1e20, -1e20, -1e20);
    for (size_t io = 0; io < OBBs.size(); io++) {
      for (size_t p = 0; p < OBBs[io].points.size(); p++) {
        // size_t i = vertexID[id];
        vec3r p_prime(r * OBBs[io].points[p], u * OBBs[io].points[p], f * OBBs[io].points[p]);
        if (minim.x > p_prime.x)
          minim.x = p_prime.x;
        if (minim.y > p_prime.y)
          minim.y = p_prime.y;
        if (minim.z > p_prime.z)
          minim.z = p_prime.z;
        if (maxim.x < p_prime.x)
          maxim.x = p_prime.x;
        if (maxim.y < p_prime.y)
          maxim.y = p_prime.y;
        if (maxim.z < p_prime.z)
          maxim.z = p_prime.z;
      }
    }

    // set the center of the OBB to be the average of the
    // minimum and maximum, and the extents be half of the
    // difference between the minimum and maximum
    fittedObb.center = eigvec * (0.5 * (maxim + minim));
    fittedObb.e[0] = r;
    fittedObb.e[1] = u;
    fittedObb.e[2] = f;
    fittedObb.extent = 0.5 * (maxim - minim);

    fittedObb.enlarge(radius); // Add the Minskowski radius

    return fittedObb;
  }

  // never call the following recursive function with empty OBBs
  // Usage:
  // OBBtree obbtree;
  // obbtree.root = OBBtree::recursiveBuild(obbtree.root, OBBs);
  static OBBnode *recursiveBuild(OBBnode *node, std::vector<OBBbundle> &OBBs, double radius = 0.0) {
    node = new OBBnode();
    node->boundary = OBBtree::fitOBB(OBBs);
    if (OBBs.size() == 1) {
      node->refId = OBBs[0].Id;
      return node;
    }

    std::vector<OBBbundle> first_subOBBs;
    std::vector<OBBbundle> second_subOBBs;
    OBBtree::subdivide(OBBs, first_subOBBs, second_subOBBs);

    if (!first_subOBBs.empty()) {
      node->first = OBBtree::recursiveBuild(node->first, first_subOBBs);
    }
    if (!second_subOBBs.empty()) {
      node->second = OBBtree::recursiveBuild(node->second, second_subOBBs);
    }
    return node;
  }

  // Usage:
  // std::vector<std::pair<size_t, size_t>> intersections;
  // OBBtree::getIntersectionsIds(obbtree1.root, obbtree2.root, intersections, posB_relativeTo_posA, QB_relativeTo_QA);
  static void getIntersectionsIds(OBBnode *nodeA, OBBnode *nodeB, std::vector<std::pair<size_t, size_t>> &intersections,
                                  double scaleFactor, vec3r &posB_relativeTo_posA, quat &QB_relativeTo_QA) {
    if (nodeA == nullptr || nodeB == nullptr) {
      return;
    }

    OBB BoundaryA = nodeA->boundary;
    BoundaryA.center *= scaleFactor;
    BoundaryA.extent *= scaleFactor;

    OBB movedBoundaryB = nodeB->boundary;
    movedBoundaryB.center *= scaleFactor;
    movedBoundaryB.extent *= scaleFactor;
    movedBoundaryB.rotate(QB_relativeTo_QA);
    movedBoundaryB.translate(posB_relativeTo_posA);

    if (BoundaryA.intersect(movedBoundaryB) == false) {
      return;
    }

    if (nodeA->isLeaf() && nodeB->isLeaf()) {
      intersections.push_back(std::pair<size_t, size_t>(nodeA->refId, nodeB->refId));
    } else if (!nodeA->isLeaf() && !nodeB->isLeaf()) {
      getIntersectionsIds(nodeA->first, nodeB->first, intersections, scaleFactor, posB_relativeTo_posA,
                          QB_relativeTo_QA);
      getIntersectionsIds(nodeA->first, nodeB->second, intersections, scaleFactor, posB_relativeTo_posA,
                          QB_relativeTo_QA);
      getIntersectionsIds(nodeA->second, nodeB->first, intersections, scaleFactor, posB_relativeTo_posA,
                          QB_relativeTo_QA);
      getIntersectionsIds(nodeA->second, nodeB->second, intersections, scaleFactor, posB_relativeTo_posA,
                          QB_relativeTo_QA);
    } else if (nodeA->isLeaf() && !nodeB->isLeaf()) {
      getIntersectionsIds(nodeA, nodeB->first, intersections, scaleFactor, posB_relativeTo_posA, QB_relativeTo_QA);
      getIntersectionsIds(nodeA, nodeB->second, intersections, scaleFactor, posB_relativeTo_posA, QB_relativeTo_QA);
    } else if (!nodeA->isLeaf() && nodeB->isLeaf()) {
      getIntersectionsIds(nodeA->first, nodeB, intersections, scaleFactor, posB_relativeTo_posA, QB_relativeTo_QA);
      getIntersectionsIds(nodeA->second, nodeB, intersections, scaleFactor, posB_relativeTo_posA, QB_relativeTo_QA);
    }
  }
};

#endif /* end of include guard: OBB_TREE_HPP */
//===========================================================================

#if 0

#include <cstdlib>
#include <iostream>

// add fake points for devel
void add6Points(OBBbundle &leaf) {
  leaf.points.push_back(leaf.obb.center - leaf.obb.extent.x * leaf.obb.e[0]);
  leaf.points.push_back(leaf.obb.center + leaf.obb.extent.x * leaf.obb.e[0]);
  leaf.points.push_back(leaf.obb.center - leaf.obb.extent.y * leaf.obb.e[1]);
  leaf.points.push_back(leaf.obb.center + leaf.obb.extent.y * leaf.obb.e[1]);
  leaf.points.push_back(leaf.obb.center - leaf.obb.extent.z * leaf.obb.e[2]);
  leaf.points.push_back(leaf.obb.center + leaf.obb.extent.z * leaf.obb.e[2]);
}

int main(int argc, char const *argv[]) {

  std::cout << "Test of OBB tree\n";

  std::vector<OBBbundle> OBBs1;
  std::vector<OBBbundle> OBBs2;

  OBB obb;

  // -----
  obb.center.set(0, 0, 0);
  obb.extent.set(1, 1, 1);
  OBBs1.push_back(OBBbundle(obb, 0));
  add6Points(OBBs1.back());

  obb.center.set(2, 2, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs1.push_back(OBBbundle(obb, 1));
  add6Points(OBBs1.back());

  obb.center.set(3, 3, 0);
  obb.extent.set(1, 1, 1);
  OBBs1.push_back(OBBbundle(obb, 2));
  add6Points(OBBs1.back());

  obb.center.set(1, 4, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs1.push_back(OBBbundle(obb, 3));
  add6Points(OBBs1.back());

  // -----

  obb.center.set(3, 4, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs2.push_back(OBBbundle(obb, 0));
  add6Points(OBBs2.back());

  obb.center.set(5, 2, 0);
  // obb.center.set(1, 1, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs2.push_back(OBBbundle(obb, 1));
  add6Points(OBBs2.back());
  // -----

  OBBtree obbtree1;
  obbtree1.root = OBBtree::recursiveBuild(obbtree1.root, OBBs1, 0.0);
  std::cout << "obbtree1.root->boundary = " << obbtree1.root->boundary << '\n';

  OBBtree obbtree2;
  obbtree2.root = OBBtree::recursiveBuild(obbtree2.root, OBBs2, 0.0);
  std::cout << "obbtree2.root->boundary = " << obbtree2.root->boundary << '\n';

  std::vector<std::pair<size_t, size_t>> intersections;
  vec3r L(-1, 0, 0);
  quat Q;
  OBBtree::getIntersectionsIds(obbtree1.root, obbtree2.root, intersections, 1.0, L, Q);
  std::cout << "intersections.size() = " << intersections.size() << '\n';

  for (size_t i = 0; i < intersections.size(); i++) {
    std::cout << intersections[i].first << " <-> " << intersections[i].second << '\n';
  }

  return 0;
}

#endif

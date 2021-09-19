#include "../OBBtree.hpp"

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
  OBBtree::TreeIntersectionIds(obbtree1.root, obbtree2.root, intersections, 1.0, L, Q);
  std::cout << "intersections.size() = " << intersections.size() << '\n';

  for (size_t i = 0; i < intersections.size(); i++) {
    std::cout << intersections[i].first << " <-> " << intersections[i].second << '\n';
  }

  return 0;
}

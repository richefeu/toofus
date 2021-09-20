#include "../OBBtree.hpp"

#include <cstdlib>
#include <iostream>

struct myData {
  int fake_data1;
  int fake_data2;
};

// add fake points for devel
void add6Points(OBBbundle<myData> &leaf) {
  leaf.points.push_back(leaf.obb.center - leaf.obb.extent.x * leaf.obb.e[0]);
  leaf.points.push_back(leaf.obb.center + leaf.obb.extent.x * leaf.obb.e[0]);
  leaf.points.push_back(leaf.obb.center - leaf.obb.extent.y * leaf.obb.e[1]);
  leaf.points.push_back(leaf.obb.center + leaf.obb.extent.y * leaf.obb.e[1]);
  leaf.points.push_back(leaf.obb.center - leaf.obb.extent.z * leaf.obb.e[2]);
  leaf.points.push_back(leaf.obb.center + leaf.obb.extent.z * leaf.obb.e[2]);
}

int main(int argc, char const *argv[]) {

  std::cout << "Test of OBB tree\n";

  std::vector<OBBbundle<myData>> OBBs1;
  std::vector<OBBbundle<myData>> OBBs2;

  OBB obb;

  // -----
  obb.center.set(0, 0, 0);
  obb.extent.set(1, 1, 1);
  OBBs1.push_back(OBBbundle<myData>(obb, {0, 28}));
  add6Points(OBBs1.back());

  obb.center.set(2, 2, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs1.push_back(OBBbundle<myData>(obb, {1, 93}));
  add6Points(OBBs1.back());

  obb.center.set(3, 3, 0);
  obb.extent.set(1, 1, 1);
  OBBs1.push_back(OBBbundle<myData>(obb, {2, 29}));
  add6Points(OBBs1.back());

  obb.center.set(1, 4, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs1.push_back(OBBbundle<myData>(obb, {3, 35}));
  add6Points(OBBs1.back());

  // -----

  OBBtree<myData> obbtree1;
  obbtree1.root = OBBtree<myData>::recursiveBuild(obbtree1.root, OBBs1, 0.0);
  std::cout << "obbtree1.root->boundary = " << obbtree1.root->boundary << '\n';

  obb.center.set(2, 5, 0);
  obb.extent.set(1.5, 1.5, 1);
  std::vector<myData> intersections;
  OBBtree<myData>::OBBIntersectionIds(obb, obbtree1.root, intersections);
  for (size_t i = 0; i < intersections.size(); i++) {
    std::cout << " -> obb sees " << intersections[i].fake_data1 << " (" << intersections[i].fake_data2 << ")\n";
  }

  return 0;
}

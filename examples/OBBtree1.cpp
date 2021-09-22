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
  OBBs1.push_back(OBBbundle<myData>(obb, {0, 56}));
  add6Points(OBBs1.back());

  obb.center.set(2, 2, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs1.push_back(OBBbundle<myData>(obb, {1, 6}));
  add6Points(OBBs1.back());

  obb.center.set(3, 3, 0);
  obb.extent.set(1, 1, 1);
  OBBs1.push_back(OBBbundle<myData>(obb, {2, 5}));
  add6Points(OBBs1.back());

  obb.center.set(1, 4, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs1.push_back(OBBbundle<myData>(obb, {3, 74}));
  add6Points(OBBs1.back());

  // -----

  obb.center.set(3, 4, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs2.push_back(OBBbundle<myData>(obb, {0, 0}));
  add6Points(OBBs2.back());

  obb.center.set(5, 2, 0);
  obb.extent.set(0.5, 0.5, 1);
  OBBs2.push_back(OBBbundle<myData>(obb, {1, 1}));
  add6Points(OBBs2.back());

  // -----

  OBBtree<myData> obbtree1;
  obbtree1.root = OBBtree<myData>::recursiveBuild(obbtree1.root, OBBs1, 0.0);
  std::cout << "obbtree1.root->boundary = " << obbtree1.root->boundary << '\n';

  OBBtree<myData> obbtree2;
  obbtree2.root = OBBtree<myData>::recursiveBuild(obbtree2.root, OBBs2, 0.0);
  std::cout << "obbtree2.root->boundary = " << obbtree2.root->boundary << '\n';

  std::vector<std::pair<myData, myData>> intersections;
  vec3r L(-1, 0, 0);
  quat Q;
  OBBtree<myData>::TreeIntersectionIds(obbtree1.root, obbtree2.root, intersections, 1.0, 1.0, 0.0, L, Q);
  std::cout << "intersections.size() = " << intersections.size() << '\n';

  for (size_t i = 0; i < intersections.size(); i++) {
    std::cout << intersections[i].first.fake_data1 << " <-> " << intersections[i].second.fake_data1 << '\n';
    std::cout << intersections[i].first.fake_data2 << "  |  " << intersections[i].second.fake_data2 << "\n\n";
  }

  return 0;
}

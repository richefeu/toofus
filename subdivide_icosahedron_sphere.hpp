#ifndef SUBDIVIDE_ICOSAHEDRON_SPHERE_HPP
#define SUBDIVIDE_ICOSAHEDRON_SPHERE_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "vec3.hpp"

template <typename T>
T subdivide_icosahedron_sphere(std::vector<vec3<T>> &points, T R, int num_subdivisions = 0,
                               std::vector<std::vector<int>> *faces = nullptr) {

#define X .525731112119133606
#define Z .850650808352039932
  const double icosahedron_vertices[12][3] = {{-X, 0.0, Z}, {X, 0.0, Z},  {-X, 0.0, -Z}, {X, 0.0, -Z},
                                              {0.0, Z, X},  {0.0, Z, -X}, {0.0, -Z, X},  {0.0, -Z, -X},
                                              {Z, X, 0.0},  {-Z, X, 0.0}, {Z, -X, 0.0},  {-Z, -X, 0.0}};
#undef X
#undef Z

  std::vector<std::vector<int>> icosahedron_faces = {{0, 4, 1},  {0, 9, 4},  {9, 5, 4},  {4, 5, 8},  {4, 8, 1},
                                                     {8, 10, 1}, {8, 3, 10}, {5, 3, 8},  {5, 2, 3},  {2, 7, 3},
                                                     {7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6}, {0, 1, 6},
                                                     {6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5},  {7, 2, 11}};

  points.clear();
  for (int i = 0; i < 12; i++) {
    vec3<T> p;
    p.x = icosahedron_vertices[i][0];
    p.y = icosahedron_vertices[i][1];
    p.z = icosahedron_vertices[i][2];
    p.normalize();
    p *= R;
    points.push_back(p);
  }

  // Subdivise chaque triangle en quatre triangles plus petits
  for (int i = 0; i < num_subdivisions; i++) {
    std::vector<std::vector<int>> new_faces;
    size_t faces_initial_size  = icosahedron_faces.size();
    size_t points_initial_size = points.size();
    for (size_t j = 0; j < faces_initial_size; j++) {
      int v1 = icosahedron_faces[j][0];
      int v2 = icosahedron_faces[j][1];
      int v3 = icosahedron_faces[j][2];

      vec3<T> p1 = points[v1];
      vec3<T> p2 = points[v2];
      vec3<T> p3 = points[v3];

      vec3<T> mid1 = 0.5 * (p1 + p2);
      vec3<T> mid2 = 0.5 * (p2 + p3);
      vec3<T> mid3 = 0.5 * (p1 + p3);

      int m1 = points.size();
      int m2 = m1 + 1;
      int m3 = m2 + 1;
      points.push_back(mid1);
      points.push_back(mid2);
      points.push_back(mid3);

      new_faces.push_back({v1, m1, m3});
      new_faces.push_back({v2, m2, m1});
      new_faces.push_back({v3, m3, m2});
      new_faces.push_back({m1, m2, m3});
    }

    // Remplacer par les nouveaux triangles
    icosahedron_faces = new_faces;

    // Normalise les points et les déplace sur la sphère de rayon R
    for (size_t i = points_initial_size; i < points.size(); i++) {
      points[i].normalize();
      points[i] *= R;
    }
  }

  if (faces != nullptr) { *faces = icosahedron_faces; }

  return norm(points[icosahedron_faces[0][0]] - points[icosahedron_faces[0][1]]);
}

#endif

#if 0

#include <iomanip>
int main() {
  double R = 0.01;
  std::vector<vec3r> points;
  std::vector<std::vector<int>> faces;
  double dst = subdivide_icosahedron_sphere(points, R, 0, &faces);

  std::ofstream file("points.txt");
  file << std::setprecision(15);
  for (const auto& p : points) {
    file << p << ' ' << dst * 0.5 << std::endl;
  }
  std::cout << "nombre de points = " << points.size() << '\n';
  std::cout << "distance des points = " << dst << '\n';
  std::cout << "nombre de face = " << faces.size() << '\n';

  return 0;
}

#endif
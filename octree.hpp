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

#ifndef OCTREE_HPP
#define OCTREE_HPP

#include <cmath>
#include <vector>

struct ot_Point {
  double x;
  double y;
  double z;
  void *userData;

  ot_Point() : x(0.0), y(0.0), z(0.0), userData(nullptr) {}
  ot_Point(double x_, double y_, double z_, void *userData_) : x(x_), y(y_), z(z_), userData(userData_) {}
};

struct ot_Sphere {
  double x;
  double y;
  double z;
  double r;
	//void *userData; // In the case where the sphere is used as a contained element

  ot_Sphere() : x(0.0), y(0.0), z(0.0), r(0.0) {}
  ot_Sphere(double x_, double y_, double z_, double r_) : x(x_), y(y_), z(z_), r(r_) {}

  bool contains(ot_Point &p) {
    double dx = p.x - x;
    double dy = p.y - y;
    double dz = p.z - z;
    double dd = dx * dx + dy * dy + dz * dz;
    double rr = r * r;
    return !(dd > rr);
  }
};

// An Axis Aligned (Bounding or not) Box
struct ot_Box {
  double xmin;
  double ymin;
  double zmin;
  double xmax;
  double ymax;
  double zmax;
	//void *userData; // In the case where the Axis Aligned Box is used as a contained element

  ot_Box() : xmin(0.0), ymin(0.0), zmin(0.0), xmax(0.0), ymax(0.0), zmax(0.0) {}
  ot_Box(double xmin_, double ymin_, double zmin_, double xmax_, double ymax_, double zmax_)
      : xmin(xmin_), ymin(ymin_), zmin(zmin_), xmax(xmax_), ymax(ymax_), zmax(zmax_) {}

  bool contains(ot_Point &p) {
    return !(p.x < xmin || p.x > xmax || p.y < ymin || p.y > ymax || p.z < zmin || p.z > zmax);
  }

  bool intersects(ot_Box &range) {
    return !(xmin > range.xmax || xmax < range.xmin || ymin > range.ymax || ymax < range.ymin || zmin > range.zmax ||
             zmax < range.zmin);
  }

  bool intersects(ot_Sphere &c) {
    double dx = c.x - std::max(xmin, std::min(c.x, xmax));
    double dy = c.y - std::max(ymin, std::min(c.y, ymax));
    double dz = c.z - std::max(zmin, std::min(c.z, zmax));
    return (dx * dx + dy * dy + dz * dz) < (c.r * c.r);
  }
};

class OcTree {
private:
  ot_Box boundary;
  size_t capacity;
  std::vector<ot_Point> points;

  bool divided;

  OcTree *xmin_ymin_zmin;
  OcTree *xmax_ymin_zmin;
  OcTree *xmax_ymax_zmin;
  OcTree *xmin_ymax_zmin;
  OcTree *xmin_ymin_zmax;
  OcTree *xmax_ymin_zmax;
  OcTree *xmax_ymax_zmax;
  OcTree *xmin_ymax_zmax;

  void subdivide() {
    double xmid = 0.5 * (boundary.xmin + boundary.xmax);
    double ymid = 0.5 * (boundary.ymin + boundary.ymax);
    double zmid = 0.5 * (boundary.zmin + boundary.zmax);
    double xmin = boundary.xmin;
    double ymin = boundary.ymin;
    double zmin = boundary.zmin;
    double xmax = boundary.xmax;
    double ymax = boundary.ymax;
    double zmax = boundary.zmax;

    xmin_ymin_zmin = new OcTree(xmin, ymin, zmin, xmid, ymid, zmid, capacity);
    xmax_ymin_zmin = new OcTree(xmid, ymin, zmin, xmax, ymid, zmid, capacity);
    xmax_ymax_zmin = new OcTree(xmid, ymid, zmin, xmax, ymax, zmid, capacity);
    xmin_ymax_zmin = new OcTree(xmin, ymid, zmin, xmid, ymax, zmid, capacity);
    xmin_ymin_zmax = new OcTree(xmin, ymin, zmid, xmid, ymid, zmax, capacity);
    xmax_ymin_zmax = new OcTree(xmid, ymin, zmid, xmax, ymid, zmax, capacity);
    xmax_ymax_zmax = new OcTree(xmid, ymid, zmid, xmax, ymax, zmax, capacity);
    xmin_ymax_zmax = new OcTree(xmin, ymid, zmid, xmid, ymax, zmax, capacity);
    divided = true;
  }

public:
  OcTree() {
    boundary.xmin = 0.0;
    boundary.ymin = 0.0;
    boundary.zmin = 0.0;
    boundary.xmax = 0.0;
    boundary.ymax = 0.0;
    boundary.zmax = 0.0;
    capacity = 0;
    divided = false;
  }

  OcTree(ot_Box &boundary_, size_t capacity_) {
    boundary = boundary_;
    capacity = capacity_;
    divided = false;
  }

  OcTree(double xmin_, double ymin_, double zmin_, double xmax_, double ymax_, double zmax_, size_t capacity_) {
    boundary.xmin = xmin_;
    boundary.ymin = ymin_;
    boundary.zmin = zmin_;
    boundary.xmax = xmax_;
    boundary.ymax = ymax_;
    boundary.zmax = zmax_;
    capacity = capacity_;
    divided = false;
  }

  bool insert(ot_Point &point) {
    if (!boundary.contains(point)) {
      return false;
    }

    if (points.size() < capacity) {
      points.push_back(point);
      return true;
    } else {
      if (!divided) {
        subdivide();
      }

      if (xmin_ymin_zmin->insert(point)) {
        return true;
      } else if (xmax_ymin_zmin->insert(point)) {
        return true;
      } else if (xmax_ymax_zmin->insert(point)) {
        return true;
      } else if (xmin_ymax_zmin->insert(point)) {
        return true;
      } else if (xmin_ymin_zmax->insert(point)) {
        return true;
      } else if (xmax_ymin_zmax->insert(point)) {
        return true;
      } else if (xmax_ymax_zmax->insert(point)) {
        return true;
      } else if (xmin_ymax_zmax->insert(point)) {
        return true;
      }
    }
    return false;
  }

  // Querying points within a AABB
  void query(ot_Box &range, std::vector<ot_Point> &found) {
    if (!boundary.intersects(range)) {
      return;
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (range.contains(points[i])) {
        found.push_back(points[i]);
      }
    }
    if (divided) {
      xmin_ymin_zmin->query(range, found);
      xmax_ymin_zmin->query(range, found);
      xmax_ymax_zmin->query(range, found);
      xmin_ymax_zmin->query(range, found);
      xmin_ymin_zmax->query(range, found);
      xmax_ymin_zmax->query(range, found);
      xmax_ymax_zmax->query(range, found);
      xmin_ymax_zmax->query(range, found);
    }
    return;
  }

  // Querying points within a sphere
  void query(ot_Sphere &crange, std::vector<ot_Point> &found) {
    if (!boundary.intersects(crange)) {
      return;
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (crange.contains(points[i])) {
        found.push_back(points[i]);
      }
    }
    if (divided) {
      xmin_ymin_zmin->query(crange, found);
      xmax_ymin_zmin->query(crange, found);
      xmax_ymax_zmin->query(crange, found);
      xmin_ymax_zmin->query(crange, found);
      xmin_ymin_zmax->query(crange, found);
      xmax_ymin_zmax->query(crange, found);
      xmax_ymax_zmax->query(crange, found);
      xmin_ymax_zmax->query(crange, found);
    }
    return;
  }


  // Periodically querying points within a sphere
  void query_periodic(ot_Sphere &crange, std::vector<ot_Point> &found) {
    double cpyX = 0.0;
    double cpyY = 0.0;
    double cpyZ = 0.0;
    bool cut = false;
    ot_Box range(crange.x - crange.r, crange.y - crange.r, crange.z - crange.r, crange.x + crange.r,
                 crange.y + crange.r, crange.z + crange.r);
    if (range.xmax > boundary.xmax && range.xmin < boundary.xmax) {
      cpyX = -1.0;
      cut = true;
    } else if (range.xmin < boundary.xmin && range.xmax > boundary.xmin) {
      cpyX = 1.0;
      cut = true;
    }
    if (range.ymax > boundary.ymax && range.ymin < boundary.ymax) {
      cpyY = -1.0;
      cut = true;
    } else if (range.ymin < boundary.ymin && range.ymax > boundary.ymin) {
      cpyY = 1.0;
      cut = true;
    }
    if (range.zmax > boundary.zmax && range.zmin < boundary.zmax) {
      cpyZ = -1.0;
      cut = true;
    } else if (range.zmin < boundary.zmin && range.zmax > boundary.zmin) {
      cpyZ = 1.0;
      cut = true;
    }

    query(crange, found);

    if (cut == true) {
      double LX = fabs(boundary.xmax - boundary.xmin);
      double LY = fabs(boundary.ymax - boundary.ymin);
      double LZ = fabs(boundary.zmax - boundary.zmin);

      int nbZero = 0;
      if (cpyX == 0.0)
        nbZero++;
      if (cpyY == 0.0)
        nbZero++;
      if (cpyZ == 0.0)
        nbZero++;

      // one single copy (intersects a face)
      if (nbZero == 2) {
        if (cpyX != 0.0) {
          ot_Sphere rg(crange.x + cpyX * LX, crange.y, crange.z, crange.r);
          query(rg, found);
        } else if (cpyY != 0.0) {
          ot_Sphere rg(crange.x, crange.y + cpyY * LY, crange.z, crange.r);
          query(rg, found);
        } else {
          ot_Sphere rg(crange.x, crange.y, crange.z + cpyZ * LZ, crange.r);
          query(rg, found);
        }
      } else if (nbZero == 1) { // 3 copies (intersects an edge)
        if (cpyZ == 0.0) {
          ot_Sphere rg1(crange.x + cpyX * LX, crange.y, crange.z, crange.r);
          query(rg1, found);
          ot_Sphere rg2(crange.x + cpyX * LX, crange.y + cpyY * LY, crange.z, crange.r);
          query(rg2, found);
          ot_Sphere rg3(crange.x, crange.y + cpyY * LY, crange.z, crange.r);
          query(rg3, found);
        } else if (cpyY == 0.0) {
          ot_Sphere rg1(crange.x + cpyX * LX, crange.y, crange.z, crange.r);
          query(rg1, found);
          ot_Sphere rg2(crange.x + cpyX * LX, crange.y, crange.z + cpyZ * LZ, crange.r);
          query(rg2, found);
          ot_Sphere rg3(crange.x, crange.y, crange.z + cpyZ * LZ, crange.r);
          query(rg3, found);
        } else {
          ot_Sphere rg1(crange.x, crange.y + cpyY * LY, crange.z, crange.r);
          query(rg1, found);
          ot_Sphere rg2(crange.x, crange.y + cpyY * LY, crange.z + cpyZ * LZ, crange.r);
          query(rg2, found);
          ot_Sphere rg3(crange.x, crange.y, crange.z + cpyZ * LZ, crange.r);
          query(rg3, found);
        }
      } else { // 7 copies (intersects a corner)
        ot_Sphere rg1(crange.x + cpyX * LX, crange.y, crange.z, crange.r);
        query(rg1, found);
        ot_Sphere rg2(crange.x + cpyX * LX, crange.y + cpyY * LY, crange.z, crange.r);
        query(rg2, found);
        ot_Sphere rg3(crange.x, crange.y + cpyY * LY, crange.z, crange.r);
        query(rg3, found);
        ot_Sphere rg4(crange.x, crange.y, crange.z + cpyZ * LZ, crange.r);
        query(rg4, found);
        ot_Sphere rg5(crange.x + cpyX * LX, crange.y, crange.z + cpyZ * LZ, crange.r);
        query(rg5, found);
        ot_Sphere rg6(crange.x + cpyX * LX, crange.y + cpyY * LY, crange.z + cpyZ * LZ, crange.r);
        query(rg6, found);
        ot_Sphere rg7(crange.x, crange.y + cpyY * LY, crange.z + cpyZ * LZ, crange.r);
        query(rg7, found);
      }
    }
  }
};

#endif /* end of include guard: OCTREE_HPP */

#if 0
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "ExecChrono.hpp"

int main(int argc, char const *argv[]) {

  struct Particle {
    double x, y, z;
    size_t id;
  };
  std::vector<Particle> particles;
  std::ofstream all("all.txt");
  for (size_t i = 0; i < 25000; i++) {
    Particle P;
    P.x = 1. + rand() / (double)RAND_MAX * 98.;
    P.y = 1. + rand() / (double)RAND_MAX * 98.;
    P.z = 1. + rand() / (double)RAND_MAX * 98.;
    P.id = i;
    particles.push_back(P);
    all << P.x << ' ' << P.y << ' ' << P.z << '\n';
  }

  OcTree q(0, 0, 0, 100, 100, 100, 4);

  for (size_t i = 0; i < particles.size(); i++) {
    ot_Point p(particles[i].x, particles[i].y, particles[i].z, &particles[i]);
    q.insert(p);
  }

  // ot_Box range(40, 40, 60, 60);
  // ot_Sphere crange(particles[100].x, particles[100].y, particles[100].z, 20);
  ot_Sphere crange(0, 0, 0, 20);
  std::vector<ot_Point> found;
  q.query_periodic(crange, found);
  // q.query_periodic(range, found);

  std::cout << "nb found = " << found.size() << "\n";
  std::ofstream ff("found.txt");
  for (size_t i = 0; i < found.size(); i++) {
    Particle *P = static_cast<Particle *>(found[i].userData);
    if (P->id >= 100)
      ff << found[i].x << ' ' << found[i].y << ' ' << found[i].z << '\n';
  }

  size_t count;

  ExecChrono C1;
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    for (size_t j = i + 1; j < particles.size(); j++) {
      if (fabs(particles[j].x - particles[i].x) < 2 && fabs(particles[j].y - particles[i].y) < 2 &&
          fabs(particles[j].z - particles[i].z) < 2)
        count++;
    }
  }
  C1.stop();
  std::cout << "count = " << count << '\n';

  ExecChrono C2;
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    ot_Box crange(particles[i].x - 2, particles[i].y - 2, particles[i].z - 2, particles[i].x + 2, particles[i].y + 2,
                  particles[i].z + 2);
    std::vector<ot_Point> found;
    q.query(crange, found);
    for (size_t f = 0; f < found.size(); f++) {
      Particle *Pj = static_cast<Particle *>(found[f].userData);
      if (Pj->id <= i)
        continue;
      if (fabs(Pj->x - particles[i].x) < 2 && fabs(Pj->y - particles[i].y) < 2 && fabs(Pj->z - particles[i].z) < 2)
        count++;
    }
  }
  C2.stop();
  std::cout << "count = " << count << '\n';

  return 0;
}

#endif

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
  size_t index;
  void *userDataPtr;

  ot_Point() : x(0.0), y(0.0), z(0.0), index(0), userDataPtr(nullptr) {}
  ot_Point(double x_, double y_, double z_, size_t index_, void *userDataPtr_)
      : x(x_), y(y_), z(z_), index(index_), userDataPtr(userDataPtr_) {}
};

struct ot_Shape {
  virtual ~ot_Shape() {}
  virtual bool contains(ot_Point &p) = 0;
  virtual void translate(double Tx, double Ty, double Tz) = 0;
  virtual double Xmin() = 0;
  virtual double Xmax() = 0;
  virtual double Ymin() = 0;
  virtual double Ymax() = 0;
  virtual double Zmin() = 0;
  virtual double Zmax() = 0;
};

struct ot_Sphere : public ot_Shape {
  double x;
  double y;
  double z;
  double r;

  ot_Sphere() : x(0.0), y(0.0), z(0.0), r(0.0) {}
  ot_Sphere(double x_, double y_, double z_, double r_) : x(x_), y(y_), z(z_), r(r_) {}

  void translate(double Tx, double Ty, double Tz) {
    x += Tx;
    y += Ty;
    z += Tz;
  }

  double Xmin() { return (x - r); }
  double Xmax() { return (x + r); }
  double Ymin() { return (y - r); }
  double Ymax() { return (y + r); }
  double Zmin() { return (z - r); }
  double Zmax() { return (z + r); }

  bool contains(ot_Point &p) {
    double dx = p.x - x;
    double dy = p.y - y;
    double dz = p.z - z;
    double dd = dx * dx + dy * dy + dz * dz;
    // double rr = r * r;
    return !(dd > r * r);
  }
};

// An Axis Aligned (Bounding or not) Box
struct ot_Box : public ot_Shape {
  double xmin;
  double ymin;
  double zmin;
  double xmax;
  double ymax;
  double zmax;

  ot_Box() : xmin(0.0), ymin(0.0), zmin(0.0), xmax(0.0), ymax(0.0), zmax(0.0) {}
  ot_Box(double xmin_, double ymin_, double zmin_, double xmax_, double ymax_, double zmax_)
      : xmin(xmin_), ymin(ymin_), zmin(zmin_), xmax(xmax_), ymax(ymax_), zmax(zmax_) {}

  ot_Box translated(double Tx, double Ty, double Tz) {
    return ot_Box(xmin + Tx, ymin + Ty, zmin + Tz, xmax + Tx, ymax + Ty, zmax + Tz);
  }

  void translate(double Tx, double Ty, double Tz) {
    xmin += Tx;
    ymin += Ty;
    zmin += Tz;
    xmax += Tx;
    ymax += Ty;
    zmax += Tz;
  }

  double Xmin() { return xmin; }
  double Xmax() { return xmax; }
  double Ymin() { return ymin; }
  double Ymax() { return ymax; }
  double Zmin() { return zmin; }
  double Zmax() { return zmax; }

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
    double xmin = boundary.xmin;
    double ymin = boundary.ymin;
    double zmin = boundary.zmin;
    double xmax = boundary.xmax;
    double ymax = boundary.ymax;
    double zmax = boundary.zmax;
    double xmid = 0.5 * (xmin + xmax);
    double ymid = 0.5 * (ymin + ymax);
    double zmid = 0.5 * (zmin + zmax);

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

  // Querying points within the boundaries (AABB) with ot_Box or ot_Sphere
  template <typename T> void query(T &range, std::vector<ot_Point> &found, size_t indexMin = 0) {
    if (!boundary.intersects(range)) {
      return;
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].index < indexMin)
        continue;
      if (range.contains(points[i])) {
        found.push_back(points[i]);
      }
    }
    if (divided) {
      xmin_ymin_zmin->query(range, found, indexMin);
      xmax_ymin_zmin->query(range, found, indexMin);
      xmax_ymax_zmin->query(range, found, indexMin);
      xmin_ymax_zmin->query(range, found, indexMin);
      xmin_ymin_zmax->query(range, found, indexMin);
      xmax_ymin_zmax->query(range, found, indexMin);
      xmax_ymax_zmax->query(range, found, indexMin);
      xmin_ymax_zmax->query(range, found, indexMin);
    }
    return;
  }

  // Periodically querying points within a box
  template <typename T> void query_periodic(T &range, std::vector<ot_Point> &found, size_t indexMin = 0) {
    double cpyX = 0.0;
    double cpyY = 0.0;
    double cpyZ = 0.0;
    bool cut = false;
    if (range.Xmax() > boundary.xmax && range.Xmin() < boundary.xmax) {
      cpyX = -1.0;
      cut = true;
    } else if (range.Xmin() < boundary.xmin && range.Xmax() > boundary.xmin) {
      cpyX = 1.0;
      cut = true;
    }
    if (range.Ymax() > boundary.ymax && range.Ymin() < boundary.ymax) {
      cpyY = -1.0;
      cut = true;
    } else if (range.Ymin() < boundary.ymin && range.Ymax() > boundary.ymin) {
      cpyY = 1.0;
      cut = true;
    }
    if (range.Zmax() > boundary.zmax && range.Zmin() < boundary.zmax) {
      cpyZ = -1.0;
      cut = true;
    } else if (range.Zmin() < boundary.zmin && range.Zmax() > boundary.zmin) {
      cpyZ = 1.0;
      cut = true;
    }

    query(range, found, indexMin);

    if (cut == true) {
      double LX = cpyX * fabs(boundary.xmax - boundary.xmin);
      double LY = cpyY * fabs(boundary.ymax - boundary.ymin);
      double LZ = cpyZ * fabs(boundary.zmax - boundary.zmin);

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
          T rg = range;
          rg.translate(LX, 0.0, 0.0);
          query(rg, found, indexMin);
        } else if (cpyY != 0.0) {
          T rg = range;
          rg.translate(0.0, LY, 0.0);
          query(rg, found, indexMin);
        } else {
          T rg = range;
          rg.translate(0.0, 0.0, LZ);
          query(rg, found, indexMin);
        }
      } else if (nbZero == 1) { // 3 copies (intersects an edge)
        if (cpyZ == 0.0) {
          T rg = range;
          rg.translate(LX, 0.0, 0.0);
          query(rg, found, indexMin);
          rg.translate(0.0, LY, 0.0);
          query(rg, found, indexMin);
          rg.translate(-LX, 0.0, 0.0);
          query(rg, found, indexMin);
        } else if (cpyY == 0.0) {
          T rg = range;
          rg.translate(LX, 0.0, 0.0);
          query(rg, found, indexMin);
          rg.translate(0.0, 0.0, LZ);
          query(rg, found, indexMin);
          rg.translate(-LX, 0.0, 0.0);
          query(rg, found, indexMin);
        } else {
          T rg = range;
          rg.translate(0.0, LY, 0.0);
          query(rg, found, indexMin);
          rg.translate(0.0, 0.0, LZ);
          query(rg, found, indexMin);
          rg.translate(0.0, -LY, 0.0);
          query(rg, found, indexMin);
        }
      } else { // 7 copies (intersects a corner)
        T rg = range;
        rg.translate(LX, 0.0, 0.0);
        query(rg, found, indexMin);
        rg.translate(0.0, LY, 0.0);
        query(rg, found, indexMin);
        rg.translate(-LX, 0.0, 0.0);
        query(rg, found, indexMin);
        rg.translate(0.0, -LY, LZ);
        query(rg, found, indexMin);
        rg.translate(LX, 0.0, 0.0);
        query(rg, found, indexMin);
        rg.translate(0.0, LY, 0.0);
        query(rg, found, indexMin);
        rg.translate(-LX, 0.0, 0.0);
        query(rg, found, indexMin);
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
  };
  std::vector<Particle> particles;
  std::ofstream all("all.txt");
  for (size_t i = 0; i < 25000; i++) {
    Particle P;
    P.x = 1.0 + rand() / (double)RAND_MAX * 98.0;
    P.y = 1.0 + rand() / (double)RAND_MAX * 98.0;
    P.z = 1.0 + rand() / (double)RAND_MAX * 98.0;
    particles.push_back(P);
    all << P.x << ' ' << P.y << ' ' << P.z << '\n';
  }

  OcTree q(0.0, 0.0, 0.0, 100.0, 100.0, 100.0, 4);

  for (size_t i = 0; i < particles.size(); i++) {
    ot_Point p(particles[i].x, particles[i].y, particles[i].z, i, &particles[i]);
    q.insert(p);
  }

  // ot_Sphere crange(particles[100].x, particles[100].y, particles[100].z, 20);
  // std::cout << particles[100].x << ' ' << particles[100].y << ' ' << particles[100].z << '\n';
  ot_Sphere crange(4, 4, 4, 20);
  //ot_Box crange(-20, -20, -20, 20, 20, 20);
  std::vector<ot_Point> found;
  q.query_periodic(crange, found, 100);

  std::cout << "nb found = " << found.size() << "\n";
  std::ofstream ff("found.txt");
  for (size_t i = 0; i < found.size(); i++) {
    ff << found[i].x << ' ' << found[i].y << ' ' << found[i].z << ' ' << found[i].index << '\n';
  }

  size_t count;

  double distSearch = 2.0;

  ExecChrono C1("Brute force, crange = box");
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    for (size_t j = i + 1; j < particles.size(); j++) {
      if (fabs(particles[j].x - particles[i].x) < distSearch && fabs(particles[j].y - particles[i].y) < distSearch &&
          fabs(particles[j].z - particles[i].z) < distSearch)
        count++;
    }
  }
  C1.stop();
  std::cout << "count = " << count << '\n';

  ExecChrono C2("Octree, crange = box");
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    ot_Box crange(particles[i].x - distSearch, particles[i].y - distSearch, particles[i].z - distSearch,
                  particles[i].x + distSearch, particles[i].y + distSearch, particles[i].z + distSearch);
    std::vector<ot_Point> found;
    q.query(crange, found, i + 1);
    count += found.size();
  }
  C2.stop();
  std::cout << "count = " << count << '\n';

  ExecChrono C3("Octree, crange = sphere");
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    ot_Sphere crange(particles[i].x, particles[i].y, particles[i].z, distSearch);
    std::vector<ot_Point> found;
    q.query(crange, found, i + 1);
    count += found.size();
  }
  C3.stop();
  std::cout << "count = " << count << '\n';

  ExecChrono C4("Octree, crange = box, periodic");
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    ot_Box crange(particles[i].x - distSearch, particles[i].y - distSearch, particles[i].z - distSearch,
                  particles[i].x + distSearch, particles[i].y + distSearch, particles[i].z + distSearch);
    std::vector<ot_Point> found;
    q.query_periodic(crange, found, i + 1);
    count += found.size();
  }
  C4.stop();
  std::cout << "count = " << count << '\n';

  ExecChrono C5("Octree, crange = sphere, periodic");
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    ot_Sphere crange(particles[i].x, particles[i].y, particles[i].z, distSearch);
    std::vector<ot_Point> found;
    q.query_periodic(crange, found, i + 1);
    count += found.size();
  }
  C5.stop();
  std::cout << "count = " << count << '\n';

  return 0;
}

#endif

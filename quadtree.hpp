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

#ifndef QUADTREE_HPP
#define QUADTREE_HPP

#include <cmath>
#include <set>
#include <vector>

struct qdt_Point {
  double x;
  double y;
  void *userData;
  size_t order;

  qdt_Point() : x(0.0), y(0.0), userData(nullptr), order(0) {}
  qdt_Point(double x_, double y_, void *userData_, size_t o) : x(x_), y(y_), userData(userData_), order(o) {}
  bool operator<(qdt_Point const &right) const { return (order < right.order); }
};

struct qdt_Circle {
  double x;
  double y;
  double r;

  qdt_Circle() : x(0.0), y(0.0), r(0.0) {}
  qdt_Circle(double x_, double y_, double r_) : x(x_), y(y_), r(r_) {}

  bool contains(const qdt_Point &p) {
    double dx = p.x - x;
    double dy = p.y - y;
    double dd = dx * dx + dy * dy;
    double rr = r * r;
    return !(dd > rr);
  }
};

struct qdt_Rectangle {
  double xmin;
  double ymin;
  double xmax;
  double ymax;

  qdt_Rectangle() : xmin(0.0), ymin(0.0), xmax(0.0), ymax(0.0) {}
  qdt_Rectangle(double xmin_, double ymin_, double xmax_, double ymax_)
      : xmin(xmin_), ymin(ymin_), xmax(xmax_), ymax(ymax_) {}

  bool contains(const qdt_Point &p) { return !(p.x < xmin || p.x > xmax || p.y < ymin || p.y > ymax); }

  bool intersects(const qdt_Rectangle &range) {
    return !(xmin > range.xmax || xmax < range.xmin || ymin > range.ymax || ymax < range.ymin);
  }

  bool intersects(const qdt_Circle &c) {
    double dx = c.x - std::max(xmin, std::min(c.x, xmax));
    double dy = c.y - std::max(ymin, std::min(c.y, ymax));
    return (dx * dx + dy * dy) < (c.r * c.r);
  }
};

class QuadTree {
private:
  qdt_Rectangle boundary;
  size_t capacity;
  std::vector<qdt_Point> points;

  bool divided;

  QuadTree *xmin_ymin;
  QuadTree *xmax_ymin;
  QuadTree *xmax_ymax;
  QuadTree *xmin_ymax;

  void subdivide() {
    double xmid = 0.5 * (boundary.xmin + boundary.xmax);
    double ymid = 0.5 * (boundary.ymin + boundary.ymax);
    double xmin = boundary.xmin;
    double ymin = boundary.ymin;
    double xmax = boundary.xmax;
    double ymax = boundary.ymax;

    xmin_ymin = new QuadTree(xmin, ymin, xmid, ymid, capacity);
    xmax_ymin = new QuadTree(xmid, ymin, xmax, ymid, capacity);
    xmax_ymax = new QuadTree(xmid, ymid, xmax, ymax, capacity);
    xmin_ymax = new QuadTree(xmin, ymid, xmid, ymax, capacity);
    divided = true;
  }

  void deleteTree(QuadTree *node) {
    if (node == nullptr)
      return;

    // first delete both subtrees
    deleteTree(node->xmin_ymin);
    deleteTree(node->xmax_ymin);
    deleteTree(node->xmax_ymax);
    deleteTree(node->xmin_ymax);

    // then delete the node
    delete node;
    node = nullptr;
  }

public:
  QuadTree() {
    boundary.xmin = 0.0;
    boundary.ymin = 0.0;
    boundary.xmax = 0.0;
    boundary.ymax = 0.0;
    capacity = 0;
    divided = false;
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
  }

  QuadTree(qdt_Rectangle &boundary_, size_t capacity_) : boundary(boundary_) {
    //boundary = boundary_;
    capacity = capacity_;
    divided = false;
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
  }

  QuadTree(double xmin_, double ymin_, double xmax_, double ymax_, size_t capacity_) {
    boundary.xmin = xmin_;
    boundary.ymin = ymin_;
    boundary.xmax = xmax_;
    boundary.ymax = ymax_;
    capacity = capacity_;
    divided = false;
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
  }

  void reset() {
    deleteTree(xmin_ymin);
    deleteTree(xmax_ymin);
    deleteTree(xmax_ymax);
    deleteTree(xmin_ymax);
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
    divided = false;
  }

  void reset(double xmin_, double ymin_, double xmax_, double ymax_, size_t capacity_) {
    boundary.xmin = xmin_;
    boundary.ymin = ymin_;
    boundary.xmax = xmax_;
    boundary.ymax = ymax_;
    capacity = capacity_;
    deleteTree(xmin_ymin);
    deleteTree(xmax_ymin);
    deleteTree(xmax_ymax);
    deleteTree(xmin_ymax);
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
    divided = false;
  }

  bool insert(qdt_Point &point) {
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

      if (xmin_ymin->insert(point)) {
        return true;
      } else if (xmax_ymin->insert(point)) {
        return true;
      } else if (xmax_ymax->insert(point)) {
        return true;
      } else if (xmin_ymax->insert(point)) {
        return true;
      }
    }
    return false;
  }

  void query(qdt_Rectangle &range, std::vector<qdt_Point> &found) {
    if (!boundary.intersects(range)) {
      return;
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (range.contains(points[i])) {
        found.push_back(points[i]);
      }
    }
    if (divided) {
      xmin_ymin->query(range, found);
      xmax_ymin->query(range, found);
      xmax_ymax->query(range, found);
      xmin_ymax->query(range, found);
    }
    return;
  }

  // to obtain a sorted vector, use this after the call of 'query':
  // std::vector found_vector(found_set.begin(), found_set.end());
  void query(qdt_Rectangle &range, std::set<qdt_Point> &found, size_t order_min = 0) {
    if (!boundary.intersects(range)) {
      return;
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].order >= order_min && range.contains(points[i])) {
        found.insert(points[i]);
      }
    }
    if (divided) {
      xmin_ymin->query(range, found, order_min);
      xmax_ymin->query(range, found, order_min);
      xmax_ymax->query(range, found, order_min);
      xmin_ymax->query(range, found, order_min);
    }
    return;
  }

  void query(qdt_Circle &crange, std::vector<qdt_Point> &found) {
    if (!boundary.intersects(crange)) {
      return;
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (crange.contains(points[i])) {
        found.push_back(points[i]);
      }
    }
    if (divided) {
      xmin_ymin->query(crange, found);
      xmax_ymin->query(crange, found);
      xmax_ymax->query(crange, found);
      xmin_ymax->query(crange, found);
    }
    return;
  }

  void query(qdt_Circle &crange, std::set<qdt_Point> &found, size_t order_min = 0) {
    if (!boundary.intersects(crange)) {
      return;
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].order >= order_min && crange.contains(points[i])) {
        found.insert(points[i]);
      }
    }
    if (divided) {
      xmin_ymin->query(crange, found, order_min);
      xmax_ymin->query(crange, found, order_min);
      xmax_ymax->query(crange, found, order_min);
      xmin_ymax->query(crange, found, order_min);
    }
    return;
  }

  void query_periodic(qdt_Rectangle &range, std::vector<qdt_Point> &found) {
    double cpyX = 0;
    double cpyY = 0;
    bool cut = false;
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

    query(range, found);

    if (cut == true) {
      double LX = fabs(boundary.xmax - boundary.xmin);
      double LY = fabs(boundary.ymax - boundary.ymin);

      // one single copy (intersects an edge)
      if (cpyX * cpyY == 0.0) {
        if (cpyX != 0.0) {
          qdt_Rectangle rg(range.xmin + cpyX * LX, range.ymin, range.xmax + cpyX * LX, range.ymax);
          query(rg, found);
        } else {
          qdt_Rectangle rg(range.xmin, range.ymin + cpyY * LY, range.xmax, range.ymax + cpyY * LY);
          query(rg, found);
        }
      } else { // 3 copies (intersects a corner)
        qdt_Rectangle rg1(range.xmin + cpyX * LX, range.ymin, range.xmax + cpyX * LX, range.ymax);
        query(rg1, found);
        qdt_Rectangle rg2(range.xmin + cpyX * LX, range.ymin + cpyY * LY, range.xmax + cpyX * LX,
                          range.ymax + cpyY * LY);
        query(rg2, found);
        qdt_Rectangle rg3(range.xmin, range.ymin + cpyY * LY, range.xmax, range.ymax + cpyY * LY);
        query(rg3, found);
      }
    }
  }

  void query_periodic(qdt_Circle &crange, std::vector<qdt_Point> &found) {
    double cpyX = 0;
    double cpyY = 0;
    bool cut = false;
    qdt_Rectangle range(crange.x - crange.r, crange.y - crange.r, crange.x + crange.r, crange.y + crange.r);
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

    query(crange, found);

    if (cut == true) {
      double LX = fabs(boundary.xmax - boundary.xmin);
      double LY = fabs(boundary.ymax - boundary.ymin);

      // one single copy
      if (cpyX * cpyY == 0.0) {
        if (cpyX != 0.0) {
          qdt_Circle rg(crange.x + cpyX * LX, crange.y, crange.r);
          query(rg, found);
        } else {
          qdt_Circle rg(crange.x, crange.y + cpyY * LY, crange.r);
          query(rg, found);
        }
      } else { // 3 copies
        qdt_Circle rg1(crange.x + cpyX * LX, crange.y, crange.r);
        query(rg1, found);
        qdt_Circle rg2(crange.x + cpyX * LX, crange.y + cpyY * LY, crange.r);
        query(rg2, found);
        qdt_Circle rg3(crange.x, crange.y + cpyY * LY, crange.r);
        query(rg3, found);
      }
    }
  }
};

#endif /* end of include guard: QUADTREE_HPP */

#if 0
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "ExecChrono.hpp"

int main(int argc, char const *argv[]) {

  struct Particle {
    double x, y;
    size_t id;
  };
  std::vector<Particle> particles;
  std::ofstream all("all.txt");
  for (size_t i = 0; i < 100000; i++) {
    Particle P;
    P.x = 1. + rand() / (double)RAND_MAX * 98.;
    P.y = 1. + rand() / (double)RAND_MAX * 98.;
    P.id = i;
    particles.push_back(P);
    all << P.x << ' ' << P.y << '\n';
  }

  QuadTree q(0, 0, 100, 100, 4);

  for (size_t i = 0; i < particles.size(); i++) {
    qdt_Point p(particles[i].x, particles[i].y, &particles[i], i);
    q.insert(p);
  }

  size_t ipart = 100;
  double dst = 3.0;
  qdt_Rectangle range(particles[ipart].x - dst, particles[ipart].y - dst, particles[ipart].x + dst,
                      particles[ipart].y + dst);
  // qdt_Circle range(particles[100].x, particles[100].y, dst);

  std::set<qdt_Point> found;
  q.query(range, found, ipart + 1);

  // q.query_periodic(range, found);

  std::cout << "nb found = " << found.size() << "\n";
  std::ofstream ff("found.txt");
  for (std::set<qdt_Point>::iterator it = found.begin(); it != found.end(); ++it) {
    Particle *P = static_cast<Particle *>(it->userData);
    size_t i = it->order;
    ff << i << ' ' << it->x << ' ' << it->y << '\n';
  }

  // with gnuplot:
  // set size square
  // plot 'all.txt' u 1:2 w p pt 6, 'found.txt' u 2:3 w p pt 7

  /*
  ExecChrono C1;
  size_t count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    for (size_t j = i + 1; j < particles.size(); j++) {
      if (fabs(particles[j].x - particles[i].x) < 2 && fabs(particles[j].y - particles[i].y) < 2)
        count++;
    }
  }
  C1.stop();
  std::cout << "count = " << count << '\n';

  ExecChrono C2;
  count = 0;
  for (size_t i = 0; i < particles.size(); i++) {
    qdt_Rectangle crange(particles[i].x - 2, particles[i].y - 2, particles[i].x + 2, particles[i].y + 2);
    std::vector<qdt_Point> found;
    q.query(crange, found);
    for (size_t f = 0; f < found.size(); f++) {
      Particle *Pj = static_cast<Particle *>(found[f].userData);
      if (Pj->id <= i)
        continue;
      if (fabs(Pj->x - particles[i].x) < 2 && fabs(Pj->y - particles[i].y) < 2)
        count++;
    }
  }
  C2.stop();
  std::cout << "count = " << count << '\n';
*/
  return 0;
}

#endif

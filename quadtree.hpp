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
  bool operator<(qdt_Point const &right) const {
    return (order < right.order);
  }
};

struct qdt_Circle {
  double x;
  double y;
  double r;

  qdt_Circle() : x(0.0), y(0.0), r(0.0) {}
  qdt_Circle(double x_, double y_, double r_) : x(x_), y(y_), r(r_) {}

  /// \brief Check if a point is inside the circle.
  ///
  /// This function checks if a point is inside the circle. The point is
  /// considered inside the circle if its distance from the center of the
  /// circle is not greater than the circle radius.
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

  /// \brief Construct a rectangle with specified bounds.
  ///
  /// This constructor initializes a rectangle using the given minimum and
  /// maximum x and y coordinates.
  ///
  /// \param xmin_ The minimum x-coordinate (left edge) of the rectangle.
  /// \param ymin_ The minimum y-coordinate (bottom edge) of the rectangle.
  /// \param xmax_ The maximum x-coordinate (right edge) of the rectangle.
  /// \param ymax_ The maximum y-coordinate (top edge) of the rectangle.
  qdt_Rectangle(double xmin_, double ymin_, double xmax_, double ymax_)
      : xmin(xmin_), ymin(ymin_), xmax(xmax_), ymax(ymax_) {}

  /// \brief Check if a point is inside the rectangle.
  ///
  /// This function determines if a given point is inside the rectangle
  /// by checking if the point's coordinates are within the rectangle's bounds.
  ///
  /// \param p The point to check.
  /// \return True if the point is inside the rectangle, false otherwise.
  bool contains(const qdt_Point &p) {
    return !(p.x < xmin || p.x > xmax || p.y < ymin || p.y > ymax);
  }

  /// \brief Check if the rectangle intersects with the given range.
  ///
  /// This function determines if the current rectangle intersects with the
  /// given range by checking if any of the rectangle's bounds are outside of
  /// the range.
  ///
  /// \param range The range to check for intersection.
  /// \return True if the rectangle intersects with the given range, false
  ///         otherwise.
  bool intersects(const qdt_Rectangle &range) {
    return !(xmin > range.xmax || xmax < range.xmin || ymin > range.ymax || ymax < range.ymin);
  }

  /// \brief Check if the rectangle intersects with the given circle.
  ///
  /// This function determines if the current rectangle intersects with the
  /// given circle by calculating the distance between the center of the circle
  /// and the closest point on the rectangle, and checking if this distance is
  /// less than the circle's radius.  This is done using the formula for the
  /// distance between a point and a line segment.
  ///
  /// \param c The circle to check for intersection.
  /// \return True if the rectangle intersects with the given circle, false
  ///         otherwise.
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

  /// \brief Subdivide the current quadtree node into four smaller quadrants.
  ///
  /// This function splits the current quadtree node into four smaller
  /// quadrants, creating four new child nodes. Each child node represents
  /// a quadrant of the current node's boundary. The division is performed
  /// by calculating the midpoint of the current boundary, and using these
  /// midpoints to define the boundaries of the new quadrants. After
  /// subdividing, the quadtree node is marked as divided.
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
    divided   = true;
  }

  /// \brief Recursively deletes all nodes in the quadtree.
  ///
  /// This function traverses the quadtree starting from the given node,
  /// recursively deleting all child nodes and then the node itself. This
  /// ensures that all dynamically allocated memory used by the quadtree
  /// is properly deallocated.
  ///
  /// \param node The root node of the subtree to delete. If nullptr, the
  ///             function does nothing.
  void deleteTree(QuadTree *node) {
    if (node == nullptr) return;

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
  /// \brief Construct a quadtree node with default values.
  ///
  /// All boundary values are set to zero, and the node is not divided.
  /// The capacity is also set to zero.
  QuadTree() {
    boundary.xmin = 0.0;
    boundary.ymin = 0.0;
    boundary.xmax = 0.0;
    boundary.ymax = 0.0;
    capacity      = 0;
    divided       = false;
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
  }

  /// \brief Construct a quadtree node with specified boundary and capacity.
  ///
  /// This constructor initializes a quadtree node using the provided
  /// rectangle boundary and capacity. The node is initially not divided,
  /// and all child nodes are set to nullptr.
  ///
  /// \param boundary_ The rectangle boundary of the quadtree node.
  /// \param capacity_ The maximum number of points the quadtree node can hold.
  QuadTree(qdt_Rectangle &boundary_, size_t capacity_) : boundary(boundary_) {
    capacity  = capacity_;
    divided   = false;
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
  }

  /// \brief Construct a quadtree node with specified boundary and capacity.
  ///
  /// This constructor initializes a quadtree node using the provided
  /// coordinates for the rectangle boundary and the maximum number of points
  /// the quadtree node can hold. The node is initially not divided,
  /// and all child nodes are set to nullptr.
  ///
  /// \param xmin_ The minimum x-coordinate of the quadtree node's boundary.
  /// \param ymin_ The minimum y-coordinate of the quadtree node's boundary.
  /// \param xmax_ The maximum x-coordinate of the quadtree node's boundary.
  /// \param ymax_ The maximum y-coordinate of the quadtree node's boundary.
  /// \param capacity_ The maximum number of points the quadtree node can hold.
  QuadTree(double xmin_, double ymin_, double xmax_, double ymax_, size_t capacity_) {
    boundary.xmin = xmin_;
    boundary.ymin = ymin_;
    boundary.xmax = xmax_;
    boundary.ymax = ymax_;
    capacity      = capacity_;
    divided       = false;
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
  }

  /// \brief Reset the quadtree node to its original state.
  ///
  /// This function deletes all child nodes, resets the divided flag to false,
  /// and sets all child pointers to nullptr. The quadtree node is then
  /// back in its original state, as if it had just been constructed.
  void reset() {
    deleteTree(xmin_ymin);
    deleteTree(xmax_ymin);
    deleteTree(xmax_ymax);
    deleteTree(xmin_ymax);
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
    divided                                       = false;
  }

  /// \brief Reset the quadtree node to its original state with new boundary and capacity.
  ///
  /// This function deletes all child nodes, resets the divided flag to false,
  /// and sets all child pointers to nullptr. The quadtree node is then
  /// back in its original state, as if it had just been constructed.
  ///
  /// \param xmin_ The minimum x-coordinate of the quadtree node's boundary.
  /// \param ymin_ The minimum y-coordinate of the quadtree node's boundary.
  /// \param xmax_ The maximum x-coordinate of the quadtree node's boundary.
  /// \param ymax_ The maximum y-coordinate of the quadtree node's boundary.
  /// \param capacity_ The maximum number of points the quadtree node can hold.
  void reset(double xmin_, double ymin_, double xmax_, double ymax_, size_t capacity_) {
    boundary.xmin = xmin_;
    boundary.ymin = ymin_;
    boundary.xmax = xmax_;
    boundary.ymax = ymax_;
    capacity      = capacity_;
    deleteTree(xmin_ymin);
    deleteTree(xmax_ymin);
    deleteTree(xmax_ymax);
    deleteTree(xmin_ymax);
    xmin_ymin = xmax_ymin = xmax_ymax = xmin_ymax = nullptr;
    divided                                       = false;
  }

  /// \brief Insert a point into the quadtree node.
  ///
  /// This function attempts to insert a point into the quadtree node. If the
  /// point is not within the node's boundary, the function returns false.
  /// Otherwise, if the node is not divided, the function divides the node and
  /// inserts the point into one of the child nodes. If the node is already
  /// divided, the function inserts the point into one of the child nodes.
  ///
  /// \param point The point to be inserted.
  ///
  /// \return true if the point was inserted, false otherwise.
  bool insert(qdt_Point &point) {
    if (!boundary.contains(point)) { return false; }

    if (points.size() < capacity) {
      points.push_back(point);
      return true;
    } else {
      if (!divided) { subdivide(); }

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

  /// \brief Find all points within a given range.
  ///
  /// This function finds all points within a given range and adds them to the
  /// found vector. If the node is divided, the function is recursively called
  /// for each child node.
  ///
  /// \param range The range to search for points in.
  /// \param found The vector of points found in the range.
  void query(qdt_Rectangle &range, std::vector<qdt_Point> &found) {
    if (!boundary.intersects(range)) { return; }
    for (size_t i = 0; i < points.size(); i++) {
      if (range.contains(points[i])) { found.push_back(points[i]); }
    }
    if (divided) {
      xmin_ymin->query(range, found);
      xmax_ymin->query(range, found);
      xmax_ymax->query(range, found);
      xmin_ymax->query(range, found);
    }
    return;
  }

  /// \brief Find all points within a given rectangle range with a minimum order.
  ///
  /// This function finds all points within a given rectangle range that have an
  /// order greater than or equal to the specified minimum order and adds them
  /// to the 'found' set. If the quadtree node is divided, the function is
  /// recursively called for each child node.
  ///
  /// \remark To obtain a sorted vector, use this after the call of 'query':
  /// \code
  /// std::vector found_vector(found_set.begin(), found_set.end());
  /// \endcode
  ///
  /// \param range The rectangle range to search for points in.
  /// \param found The set of points found in the range with the specified order.
  /// \param order_min The minimum order of points to be included in the result.
  void query(qdt_Rectangle &range, std::set<qdt_Point> &found, size_t order_min = 0) {
    if (!boundary.intersects(range)) { return; }
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].order >= order_min && range.contains(points[i])) { found.insert(points[i]); }
    }
    if (divided) {
      xmin_ymin->query(range, found, order_min);
      xmax_ymin->query(range, found, order_min);
      xmax_ymax->query(range, found, order_min);
      xmin_ymax->query(range, found, order_min);
    }
    return;
  }

  /// \brief Find all points within a given circle range.
  ///
  /// This function finds all points within a given circle range and adds them
  /// to the 'found' vector. If the quadtree node is divided, the function is
  /// recursively called for each child node.
  ///
  /// \param crange The circle range to search for points in.
  /// \param found The vector of points found in the range.
  void query(qdt_Circle &crange, std::vector<qdt_Point> &found) {
    if (!boundary.intersects(crange)) { return; }
    for (size_t i = 0; i < points.size(); i++) {
      if (crange.contains(points[i])) { found.push_back(points[i]); }
    }
    if (divided) {
      xmin_ymin->query(crange, found);
      xmax_ymin->query(crange, found);
      xmax_ymax->query(crange, found);
      xmin_ymax->query(crange, found);
    }
    return;
  }

  /// \brief Find all points within a given circle range with a minimum order.
  ///
  /// This function finds all points within a given circle range that have an
  /// order greater than or equal to the specified minimum order and adds them
  /// to the 'found' set. If the quadtree node is divided, the function is
  /// recursively called for each child node.
  ///
  /// \remark To obtain a sorted vector, use this after the call of 'query':
  /// \code
  /// std::vector found_vector(found_set.begin(), found_set.end());
  /// \endcode
  ///
  /// \param crange The circle range to search for points in.
  /// \param found The set of points found in the range with the specified order.
  /// \param order_min The minimum order of points to be included in the result.
  void query(qdt_Circle &crange, std::set<qdt_Point> &found, size_t order_min = 0) {
    if (!boundary.intersects(crange)) { return; }
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].order >= order_min && crange.contains(points[i])) { found.insert(points[i]); }
    }
    if (divided) {
      xmin_ymin->query(crange, found, order_min);
      xmax_ymin->query(crange, found, order_min);
      xmax_ymax->query(crange, found, order_min);
      xmin_ymax->query(crange, found, order_min);
    }
    return;
  }

  /// \brief Find all points within a given rectangle range in a periodic domain.
  ///
  /// This function finds all points within a given rectangle range in a periodic
  /// domain and adds them to the 'found' vector. If the node is divided, the
  /// function is recursively called for each child node.
  ///
  /// \param range The rectangle range to search for points in.
  /// \param found The vector of points found in the range.
  void query_periodic(qdt_Rectangle &range, std::vector<qdt_Point> &found) {
    double cpyX = 0;
    double cpyY = 0;
    bool cut    = false;
    if (range.xmax > boundary.xmax && range.xmin < boundary.xmax) {
      cpyX = -1.0;
      cut  = true;
    } else if (range.xmin < boundary.xmin && range.xmax > boundary.xmin) {
      cpyX = 1.0;
      cut  = true;
    }
    if (range.ymax > boundary.ymax && range.ymin < boundary.ymax) {
      cpyY = -1.0;
      cut  = true;
    } else if (range.ymin < boundary.ymin && range.ymax > boundary.ymin) {
      cpyY = 1.0;
      cut  = true;
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

  /// \brief Find all points within a given circle range in a periodic domain.
  ///
  /// This function finds all points within a given circle range in a periodic
  /// domain and adds them to the 'found' vector. If the node is divided, the
  /// function is recursively called for each child node.
  ///
  /// \param crange The circle range to search for points in.
  /// \param found The vector of points found in the range.
  void query_periodic(qdt_Circle &crange, std::vector<qdt_Point> &found) {
    double cpyX = 0;
    double cpyY = 0;
    bool cut    = false;
    qdt_Rectangle range(crange.x - crange.r, crange.y - crange.r, crange.x + crange.r, crange.y + crange.r);
    if (range.xmax > boundary.xmax && range.xmin < boundary.xmax) {
      cpyX = -1.0;
      cut  = true;
    } else if (range.xmin < boundary.xmin && range.xmax > boundary.xmin) {
      cpyX = 1.0;
      cut  = true;
    }
    if (range.ymax > boundary.ymax && range.ymin < boundary.ymax) {
      cpyY = -1.0;
      cut  = true;
    } else if (range.ymin < boundary.ymin && range.ymax > boundary.ymin) {
      cpyY = 1.0;
      cut  = true;
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

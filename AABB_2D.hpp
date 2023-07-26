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

#ifndef AABB_2D_HPP
#define AABB_2D_HPP

#include <vector>

#include "vec2.hpp"

// =======================================================
// Axis Aligned Bounding Box
// =======================================================
class AABB_2D {
public:
  vec2r min, max;

  // Constructors
  AABB_2D() : min(), max() {}
  explicit AABB_2D(const vec2r &v) : min(v), max(v) {}
  AABB_2D(const vec2r &v1, const vec2r &v2) : min(component_min(v1, v2)), max(component_max(v1, v2)) {}
  AABB_2D(const AABB_2D &aabb) : min(aabb.min), max(aabb.max) {}
  AABB_2D(double xmin, double xmax, double ymin, double ymax) {
    min.set(xmin, ymin);
    max.set(xmax, ymax);
  }

  explicit AABB_2D(const std::vector<vec2r> &cloud) : min(cloud[0]), max(cloud[0]) {
    for (size_t i = 1; i < cloud.size(); ++i) {
      min = component_min(min, cloud[i]);
      max = component_max(max, cloud[i]);
    }
  }

  AABB_2D &operator=(const AABB_2D &aabb) {
    min = aabb.min;
    max = aabb.max;
    return (*this);
  }

  /**
    @brief Get the radius of a sphere that surrounds the AABB_2D,
           centered at the AABB_2D center
  */
  double getRadius() const { return 0.25 * (max - min).length(); }

  vec2r getCenter() const { return 0.5 * (min + max); }

  /**
    @brief The AABB_2D is set to a single point
  */
  void set_single(const vec2r &v) {
    min = v;
    max = v;
  }

  void add(const vec2r &v) {
    min = component_min(min, v);
    max = component_max(max, v);
  }

  void enlarge(double more) {
    min.x -= more;
    min.y -= more;
    max.x += more;
    max.y += more;
  }

  void enlarge(const vec2r &more) {
    min.x -= more.x;
    min.y -= more.y;
    max.x += more.x;
    max.y += more.y;
  }

  void merge(const AABB_2D &more) {
    min = component_min(min, more.min);
    max = component_max(max, more.max);
  }
  void enlarge(const AABB_2D &more) { // for compability
    min = component_min(min, more.min);
    max = component_max(max, more.max);
  }

  void translate(const vec2r &v) {
    min += v;
    max += v;
  }

  bool intersect(const AABB_2D &a) const {
    if (max.x < a.min.x || a.max.x < min.x || max.y < a.min.y || a.max.y < min.y)
      return false;
    return true;
  }

  bool intersect(const vec2r &a) const {
    if (max.x < a.x || a.x < min.x || max.y < a.y || a.y < min.y)
      return false;
    return true;
  }

  /**
    @brief Check intersection only in the X direction
  */
  bool intersectX(const AABB_2D &a) const {
    if (max.x < a.min.x || a.max.x < min.x)
      return false;
    return true;
  }

  /**
    @brief Check intersection only in the Y direction
  */
  bool intersectY(const AABB_2D &a) const {
    if (max.y < a.min.y || a.max.y < min.y)
      return false;
    return true;
  }
};

#endif /* end of include guard: AABB_2D_HPP */

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

#ifndef AABB_HPP
#define AABB_HPP

#include <vector>

#include "vec3.hpp"

// =======================================================
// Axis Aligned Bounding Box
// =======================================================
class AABB {
public:
  vec3r min, max;

  // Constructors
  AABB() : min(), max() {}
  explicit AABB(const vec3r &v) : min(v), max(v) {}
  AABB(const vec3r &v1, const vec3r &v2) : min(component_min(v1, v2)), max(component_max(v1, v2)) {}
  AABB(const AABB &aabb) : min(aabb.min), max(aabb.max) {}

  explicit AABB(const std::vector<vec3r> &cloud) : min(cloud[0]), max(cloud[0]) {
    for (size_t i = 1; i < cloud.size(); ++i) {
      min = component_min(min, cloud[i]);
      max = component_max(max, cloud[i]);
    }
  }

  AABB &operator=(const AABB &aabb) {
    min = aabb.min;
    max = aabb.max;
    return (*this);
  }

  /**
    @brief Get the radius of a sphere that surrounds the AABB, 
           centered at the AABB center
  */
  double getRadius() const { return 0.25 * (max - min).length(); }

  /**
    @brief The AABB is set to a single point
  */
  void set_single(const vec3r &v) {
    min = v;
    max = v;
  }

  
  void add(const vec3r &v) {
    min = component_min(min, v);
    max = component_max(max, v);
  }

  void enlarge(double more) {
    min.x -= more;
    min.y -= more;
    min.z -= more;
    max.x += more;
    max.y += more;
    max.z += more;
  }

  void enlarge(const vec3r &more) {
    min.x -= more.x;
    min.y -= more.y;
    min.z -= more.z;
    max.x += more.x;
    max.y += more.y;
    max.z += more.z;
  }

  /**
    @brief Join the 2 AABB 
    TODO : rename it 'union' or something like this
  */
  void enlarge(const AABB &more) {
    min = component_min(min, more.min);
    max = component_max(max, more.max);
  }

  void translate(const vec3r &v) {
    min += v;
    max += v;
  }

  bool intersect(const AABB &a) const {
    if (max.x < a.min.x || a.max.x < min.x || max.y < a.min.y || a.max.y < min.y || max.z < a.min.z || a.max.z < min.z)
      return false;
    return true;
  }

  bool intersect(const vec3r &a) const {
    if (max.x < a.x || a.x < min.x || max.y < a.y || a.y < min.y || max.z < a.z || a.z < min.z)
      return false;
    return true;
  }

  bool intersectX(const AABB &a) const {
    if (max.x < a.min.x || a.max.x < min.x)
      return false;
    return true;
  }

  bool intersectY(const AABB &a) const {
    if (max.y < a.min.y || a.max.y < min.y)
      return false;
    return true;
  }

  bool intersectZ(const AABB &a) const {
    if (max.z < a.min.z || a.max.z < min.z)
      return false;
    return true;
  }
};

#endif /* end of include guard: AABB_HPP */

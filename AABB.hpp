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

  /**
    @brief Constructs the AABB to encompass all points in a given cloud.

    @param points A vector of vec3r points used to determine the bounding box.
                  Initializes the min and max with the first point and iteratively
                  updates to include all points in the cloud.
   */
  explicit AABB(const std::vector<vec3r> &points) : min(points[0]), max(points[0]) {
    for (size_t i = 1; i < points.size(); ++i) {
      min = component_min(min, points[i]);
      max = component_max(max, points[i]);
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

  vec3r getCenter() const { return 0.5 * (min + max); }

  /**
    @brief The AABB is set to a single point
  */
  void set_single(const vec3r &v) {
    min = v;
    max = v;
  }

  /**
    @brief Expands the AABB to include a given point.

    @param v A vec3r point used to update the bounding box.
            Adjusts the min and max to ensure the point is within the AABB.
  */
  void add(const vec3r &v) {
    min = component_min(min, v);
    max = component_max(max, v);
  }

  /**
    @brief Expands the AABB by a given amount in all directions.

    @param more A double value used to increase the size of the AABB.
              The min coordinates are decremented by more, and the max
              coordinates are incremented by more.
  */
  void enlarge(double more) {
    min.x -= more;
    min.y -= more;
    min.z -= more;
    max.x += more;
    max.y += more;
    max.z += more;
  }

  /**
    @brief Expands the AABB by the given amounts in the x, y, and z directions.

    @param more A vec3r object used to increase the size of the AABB.
              The min coordinates are decremented by the corresponding
              values in more, and the max coordinates are incremented
              by the corresponding values in more.
  */
  void enlarge(const vec3r &more) {
    min.x -= more.x;
    min.y -= more.y;
    min.z -= more.z;
    max.x += more.x;
    max.y += more.y;
    max.z += more.z;
  }

  /**
    @brief Merges the current AABB with another AABB.

    @param more The AABB to merge with the current one. The min and max
                values of the current AABB are updated to encompass both
                AABBs.
  */
  void merge(const AABB &more) {
    min = component_min(min, more.min);
    max = component_max(max, more.max);
  }

  /**
    @brief Expands the AABB to encompass another AABB.

    @param more The AABB to merge with the current one. The min and max
              values of the current AABB are updated to encompass both
              AABBs.
  */
  void enlarge(const AABB &more) { // for compability
    min = component_min(min, more.min);
    max = component_max(max, more.max);
  }

  /**
    @brief Translates the AABB by a given vec3r vector.

    @param v The vec3r vector used to translate the AABB. The min and max
            coordinates of the AABB are incremented by the corresponding
            values in the vector.
  */
  void translate(const vec3r &v) {
    min += v;
    max += v;
  }

  /**
    @brief Checks if the current AABB intersects with another AABB.

    @param a The AABB to intersect with the current one. The AABBs are
            considered to intersect if they have a non-zero volume in
            common.

    @returns true if the AABBs intersect, false otherwise.
  */
  bool intersect(const AABB &a) const {
    if (max.x < a.min.x || a.max.x < min.x || max.y < a.min.y || a.max.y < min.y || max.z < a.min.z || a.max.z < min.z)
      return false;
    return true;
  }

  /**
    @brief Checks if the current AABB intersects with a given point.

    @param a The vec3r point to intersect with the current AABB. The AABB
            is considered to intersect with the point if the point is
            contained within the AABB.

    @returns true if the AABB intersects with the point, false otherwise.
  */
  bool intersect(const vec3r &a) const {
    if (max.x < a.x || a.x < min.x || max.y < a.y || a.y < min.y || max.z < a.z || a.z < min.z)
      return false;
    return true;
  }

  /**
    @brief Check intersection only in the X direction
  */
  bool intersectX(const AABB &a) const {
    if (max.x < a.min.x || a.max.x < min.x)
      return false;
    return true;
  }

  /**
    @brief Check intersection only in the Y direction
  */
  bool intersectY(const AABB &a) const {
    if (max.y < a.min.y || a.max.y < min.y)
      return false;
    return true;
  }

  /**
    @brief Check intersection only in the Z direction
  */
  bool intersectZ(const AABB &a) const {
    if (max.z < a.min.z || a.max.z < min.z)
      return false;
    return true;
  }
};

#endif /* end of include guard: AABB_HPP */

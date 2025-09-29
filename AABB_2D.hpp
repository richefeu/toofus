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

  /**
    @brief Constructs the AABB_2D to encompass all points in a given cloud.

    @param cloud A vector of vec2r points used to determine the bounding box.
                 Initializes the min and max with the first point and iteratively
                 updates to include all points in the cloud.
  */
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
  double getRadius() const {
    return 0.25 * (max - min).length();
  }

  vec2r getCenter() const {
    return 0.5 * (min + max);
  }

  /**
    @brief The AABB_2D is set to a single point
  */
  void set_single(const vec2r &v) {
    min = v;
    max = v;
  }

  /**
    @brief Expands the AABB_2D to include a given point.

    @param v A vec2r point used to update the bounding box.
             Adjusts the min and max to ensure the point is within the AABB_2D.
  */
  void add(const vec2r &v) {
    min = component_min(min, v);
    max = component_max(max, v);
  }

  /**
    @brief Expands the AABB_2D by a given amount in all directions.

    @param more A double value used to increase the size of the AABB_2D.
               The min coordinates are decremented by more, and the max
               coordinates are incremented by more.
  */
  void enlarge(double more) {
    min.x -= more;
    min.y -= more;
    max.x += more;
    max.y += more;
  }

  /**
    @brief Expands the AABB_2D by a specified amount in each direction.

    @param more A vec2r object specifying the amount to expand the AABB_2D.
                The min coordinates are decremented by the respective components of more,
                and the max coordinates are incremented by the respective components of more.
  */
  void enlarge(const vec2r &more) {
    min.x -= more.x;
    min.y -= more.y;
    max.x += more.x;
    max.y += more.y;
  }

  /**
    @brief Merges the current AABB_2D with another AABB_2D.

    @param more The AABB_2D to merge with the current one. The min and max
                values of the current AABB_2D are updated to encompass both
                AABB_2Ds.
  */
  void merge(const AABB_2D &more) {
    min = component_min(min, more.min);
    max = component_max(max, more.max);
  }

  /**
    @brief Expands the AABB_2D to include the given AABB_2D.

    @param more The AABB_2D to include in the current AABB_2D.
                The min and max values of the current AABB_2D are updated
                to ensure that the given AABB_2D is fully enclosed.
  */
  void enlarge(const AABB_2D &more) { // for compability
    min = component_min(min, more.min);
    max = component_max(max, more.max);
  }

  /**
    @brief Translates the AABB_2D by a given vec2r vector.

    @param v The vec2r vector used to translate the AABB_2D. The min and max
             coordinates of the AABB_2D are incremented by the corresponding
             values in the vector.
  */
  void translate(const vec2r &v) {
    min += v;
    max += v;
  }

  /**
    @brief Checks if the current AABB_2D intersects with another AABB_2D.

    @param a The AABB_2D to intersect with the current one. The AABB_2Ds are
             considered to intersect if they have a non-zero area in common.

    @returns true if the AABB_2Ds intersect, false otherwise.
  */
  bool intersect(const AABB_2D &a) const {
    if (max.x < a.min.x || a.max.x < min.x || max.y < a.min.y || a.max.y < min.y) return false;
    return true;
  }

  /**
    @brief Checks if the current AABB_2D intersects with a given point.

    @param a The point to intersect with the current AABB_2D. The AABB_2D is
             considered to intersect with the point if the point is inside the
             AABB_2D.

    @returns true if the AABB_2D intersects with the point, false otherwise.
  */
  bool intersect(const vec2r &a) const {
    if (max.x < a.x || a.x < min.x || max.y < a.y || a.y < min.y) return false;
    return true;
  }

  /**
    @brief Check intersection only in the X direction
  */
  bool intersectX(const AABB_2D &a) const {
    if (max.x < a.min.x || a.max.x < min.x) return false;
    return true;
  }

  /**
    @brief Check intersection only in the Y direction
  */
  bool intersectY(const AABB_2D &a) const {
    if (max.y < a.min.y || a.max.y < min.y) return false;
    return true;
  }
};

#endif /* end of include guard: AABB_2D_HPP */

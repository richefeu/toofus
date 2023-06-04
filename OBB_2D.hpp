#ifndef OBB_2D_HPP
#define OBB_2D_HPP

#include <cmath>
#include <limits>

#include "vec2.hpp"

struct OBB_2D {

  vec2r position;  // coordinate of center
  double width;    // Width of the OBB
  double height;   // Height of the OBB
  double rotation; // Rotation angle in radians

  OBB_2D(double t_x, double t_y, double t_width, double t_height, double t_rotation)
      : position(t_x, t_y), width(t_width), height(t_height), rotation(t_rotation) {}

  void enlarge(double more) {
    double inc = 2.0 * more;
    width += inc;
    height += inc;
  }

  bool intersect(const OBB_2D &other) const {
    // Calculate the vertices of both OBBs
    std::array<vec2r, 4> vertices1 = calculateVertices();
    std::array<vec2r, 4> vertices2 = other.calculateVertices();

    // Calculate the axes to test for separation
    std::array<vec2r, 4> axes1;
    calculateAxes(vertices1, axes1);

    // Check for separation along each axis
    for (const auto &axis : axes1) {
      if (isSeparatedByAxis(vertices1, vertices2, axis))
        return false;
    }
		
    std::array<vec2r, 4> axes2;
		calculateAxes(vertices2, axes2);
		
    // Check for separation along each axis
    for (const auto &axis : axes2) {
      if (isSeparatedByAxis(vertices1, vertices2, axis))
        return false;
    }

    // No separation along any axis, the OBBs intersect
    return true;
  }

private:
  std::array<vec2r, 4> calculateVertices() const {
    // Calculate the half extents
    double halfWidth = width / 2.0;
    double halfHeight = height / 2.0;

    // Calculate the rotation matrix
    double c = cos(rotation);
    double s = sin(rotation);

    // Calculate the local axes
    vec2r xAxis(c, s);
    vec2r yAxis(-s, c);

    // Calculate the vertices
    // vec2r center(position.x, position.y);
    std::array<vec2r, 4> vertices;
    vertices[0] = position + xAxis * halfWidth + yAxis * halfHeight;
    vertices[1] = position - xAxis * halfWidth + yAxis * halfHeight;
    vertices[2] = position - xAxis * halfWidth - yAxis * halfHeight;
    vertices[3] = position + xAxis * halfWidth - yAxis * halfHeight;

    return vertices;
  }

  void calculateAxes(const std::array<vec2r, 4> &vertices, std::array<vec2r, 4> &axes) const {

    vec2r edge;

    edge = vertices[1] - vertices[0];
    //edge.normalize();
    axes[0].set(-edge.y, edge.x);

    edge = vertices[2] - vertices[1];
    //edge.normalize();
    axes[1].set(-edge.y, edge.x);

    edge = vertices[3] - vertices[2];
    //edge.normalize();
    axes[2].set(-edge.y, edge.x);

    edge = vertices[0] - vertices[3];
    //edge.normalize();
    axes[3].set(-edge.y, edge.x);
  }

  bool isSeparatedByAxis(const std::array<vec2r, 4> &vertices1, const std::array<vec2r, 4> &vertices2,
                         const vec2r &axis) const {
    double minProjection1 = std::numeric_limits<double>::max();
    double maxProjection1 = std::numeric_limits<double>::lowest();
    double minProjection2 = std::numeric_limits<double>::max();
    double maxProjection2 = std::numeric_limits<double>::lowest();

    for (const auto &vertex : vertices1) {
      double projection = vertex * axis;
      minProjection1 = std::min(minProjection1, projection);
      maxProjection1 = std::max(maxProjection1, projection);
    }

    for (const auto &vertex : vertices2) {
      double projection = vertex * axis;
      minProjection2 = std::min(minProjection2, projection);
      maxProjection2 = std::max(maxProjection2, projection);
    }

    return !(maxProjection1 >= minProjection2 && maxProjection2 >= minProjection1);
  }
};

#endif /* end of include guard: OBB_2D_HPP */

#if 0

#include <iostream>
int main() {

  OBB_2D obb1(4.0, 4.0, 4.0, 4.0, 0.0);
  OBB_2D obb2(0.0, 0.0, 32.0, 1.0, 0.0);

  for (double a = 0.0; a < 3.14159; a += 0.04) {
    obb2.rotation = a;
    bool result = obb2.intersect(obb1);
    std::cout << "a = " << a << ", Intersect: " << (result ? "true" : "false") << std::endl;
  }

  return 0;
}

#endif

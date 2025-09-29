#ifndef OBB_2D_HPP
#define OBB_2D_HPP

#include <cmath>
#include <limits>

#include "message.hpp"
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

  bool intersects(const OBB_2D &other) const {
    // Calculate the vertices of both OBBs
    std::array<vec2r, 4> vertices1 = calculateVertices();
    std::array<vec2r, 4> vertices2 = other.calculateVertices();

    vec2r axis;

    double c = cos(rotation);
    double s = sin(rotation);
    axis.set(c, s);
    if (isSeparatedByAxis(vertices1, vertices2, axis)) { return false; }

    axis.set(-s, c);
    if (isSeparatedByAxis(vertices1, vertices2, axis)) { return false; }

    c = cos(other.rotation);
    s = sin(other.rotation);
    axis.set(c, s);
    if (isSeparatedByAxis(vertices1, vertices2, axis)) { return false; }

    axis.set(-s, c);
    if (isSeparatedByAxis(vertices1, vertices2, axis)) { return false; }

    // No separation along any axis, the OBBs intersect
    return true;
  }

  bool contains(const vec2r &point) const {
    // Calculate the half extents
    double halfWidth  = width / 2.0;
    double halfHeight = height / 2.0;

    // Calculate the rotation matrix
    double c = cos(rotation);
    double s = sin(rotation);

    // Calculate the local axes
    vec2r WAxis(c, s);
    vec2r HAxis(-s, c);

    // Calculate the point relative to the OBB center
    vec2r localPoint = point - position;

    // Calculate the projected distances onto each local axis
    double projectedX = localPoint * WAxis;
    double projectedY = localPoint * HAxis;

    // Check if the point is within the OBB's boundaries
    return (abs(projectedX) <= halfWidth) && (abs(projectedY) <= halfHeight);
  }

  bool contains(const OBB_2D &other) const {
    return contains(other.position);
  }

private:
  std::array<vec2r, 4> calculateVertices() const {
    // Calculate the half extents
    double halfWidth  = width / 2.0;
    double halfHeight = height / 2.0;

    // Calculate the rotation matrix
    double c = cos(rotation);
    double s = sin(rotation);

    // Calculate the local axes
    vec2r WAxis(c, s);
    vec2r HAxis(-s, c);

    // Calculate the vertices
    std::array<vec2r, 4> vertices;
    vertices[0] = position + WAxis * halfWidth + HAxis * halfHeight;
    vertices[1] = position - WAxis * halfWidth + HAxis * halfHeight;
    vertices[2] = position - WAxis * halfWidth - HAxis * halfHeight;
    vertices[3] = position + WAxis * halfWidth - HAxis * halfHeight;

    return vertices;
  }

  bool isSeparatedByAxis(const std::array<vec2r, 4> &vertices1, const std::array<vec2r, 4> &vertices2,
                         const vec2r &axis) const {
    double minProjection1 = 1.0e20;
    double maxProjection1 = -1.0e20;
    double minProjection2 = 1.0e20;
    double maxProjection2 = -1.0e20;

    for (const auto &vertex : vertices1) {
      double projection = vertex * axis;
      minProjection1    = std::min(minProjection1, projection);
      maxProjection1    = std::max(maxProjection1, projection);
    }

    for (const auto &vertex : vertices2) {
      double projection = vertex * axis;
      minProjection2    = std::min(minProjection2, projection);
      maxProjection2    = std::max(maxProjection2, projection);
    }

    return (maxProjection1 < minProjection2 || maxProjection2 < minProjection1);
  }
};

#endif /* end of include guard: OBB_2D_HPP */

#if 0

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Mth.hpp"
#include "message.hpp"

void generateSVG(const char *fname, OBB_2D &obbBlue, const std::vector<OBB_2D> &obbsBlack,
                 const std::vector<OBB_2D> &obbsRed, double xMin, double yMin, double xMax, double yMax) {
  double scaleFactor = 800.0 / std::max(xMax - xMin, yMax - yMin);

  std::ofstream file(fname);
  if (!file.is_open()) {
    std::cout << "Failed to open " << fname << " file" << std::endl;
    return;
  }

  // Calculate the adjusted coordinates of the OBBs
  auto adjustCoordinates = [&](double x, double y) {
    double adjustedX = (x - xMin) * scaleFactor;
    double adjustedY = ((yMax - yMin) - (y - yMin)) * scaleFactor;
    return std::make_pair(adjustedX, adjustedY);
  };

  // SVG header
  file << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"800\" height=\"800\">" << std::endl;

  file << std::setprecision(15);

  // Draw the blue OBB
  double adjustedX, adjustedY;
  std::tie(adjustedX, adjustedY) = adjustCoordinates(obbBlue.position.x, obbBlue.position.y);
  double adjustedWidth = obbBlue.width * scaleFactor;
  double adjustedHeight = obbBlue.height * scaleFactor;
  double adjustedAngle = obbBlue.rotation;
  file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2 << "\" width=\""
       << adjustedWidth << "\" height=\"" << adjustedHeight << "\" transform=\"rotate(" << -adjustedAngle * 180.0 / M_PI
       << " " << adjustedX << " " << adjustedY << ")\" fill=\"blue\" stroke=\"blue\" />" << std::endl;

  // Draw OBBs in black color
  for (const OBB_2D &obb : obbsBlack) {
    double adjustedX, adjustedY;
    std::tie(adjustedX, adjustedY) = adjustCoordinates(obb.position.x, obb.position.y);
    double adjustedWidth = obb.width * scaleFactor;
    double adjustedHeight = obb.height * scaleFactor;
    double adjustedAngle = obb.rotation;

    file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2
         << "\" width=\"" << adjustedWidth << "\" height=\"" << adjustedHeight << "\" transform=\"rotate("
         << -adjustedAngle * 180.0 / M_PI << " " << adjustedX << " " << adjustedY
         << ")\" fill=\"none\" stroke=\"black\" />" << std::endl;
  }

  // Draw OBBs in red color
  for (const OBB_2D &obb : obbsRed) {
    double adjustedX, adjustedY;
    std::tie(adjustedX, adjustedY) = adjustCoordinates(obb.position.x, obb.position.y);
    double adjustedWidth = obb.width * scaleFactor;
    double adjustedHeight = obb.height * scaleFactor;
    double adjustedAngle = obb.rotation;

    file << "<rect x=\"" << adjustedX - adjustedWidth / 2 << "\" y=\"" << adjustedY - adjustedHeight / 2
         << "\" width=\"" << adjustedWidth << "\" height=\"" << adjustedHeight << "\" transform=\"rotate("
         << -adjustedAngle * 180.0 / M_PI << " " << adjustedX << " " << adjustedY
         << ")\" fill=\"none\" stroke=\"red\" />" << std::endl;
  }

  // SVG footer
  file << "</svg>" << std::endl;

  file.close();
  std::cout << "SVG file generated: " << fname << std::endl;
}

#endif

#if 0
#include <iostream>
int main() {

  OBB_2D obb1(4.0, 4.0, 4.0, 4.0, 0.0);
  OBB_2D obb2(0.0, 0.0, 32.0, 1.0, 0.0);

  for (double a = 0.0; a < 3.14159; a += 0.04) {
    obb2.rotation = a;
    bool result = obb2.intersects(obb1);
    std::cout << "a = " << a << ", Intersect: " << (result ? "true" : "false") << std::endl;
  }

  return 0;
}
#endif

#if 0
#include <iostream>
int main() {

  OBB_2D obb1(0.0, 0.0, 2.0, 2.0, 0.2);
  OBB_2D obb2(1.8, 1.99, 2.0, 2.0, -0.1);

  bool result = obb2.intersects(obb1);
  std::cout << "Intersect: " << (result ? "true" : "false") << std::endl;

  std::vector<OBB_2D> obbsBlack, obbsRed;
  if (result)
    obbsRed.push_back(obb2);
  else
    obbsBlack.push_back(obb2);

  generateSVG("2OBBs.svg", obb1, obbsBlack, obbsRed, -3, -3, 4, 4);

  return 0;
}
#endif

#if 0
#include <iostream>
int main() {

  OBB_2D obb(0.0, 0.0, 4.0, 4.0, 0.0);
  vec2r P(2.1, 2.);
  bool result = obb.contains(P);
  std::cout << "Contains point (" << P << "): " << (result ? "true" : "false") << std::endl;

  OBB_2D obb2(1.3, -1.0, 2.0, 2.0, 0.78);
  result = obb.contains(obb2);
  std::cout << "Contains obb at position (" << obb2.position << "): " << (result ? "true" : "false") << std::endl;

  return 0;
}
#endif

#if 0

int main() {

  // Set up random number generation
  std::vector<OBB_2D> obbs;
  for (double centerX = -16; centerX <= 16; centerX += 1.0) {
    for (double centerY = -16; centerY <= 16; centerY += 1.0) {
      double width = Mth::random(0.6, 1.2);
      double height = Mth::random(0.6, 1.2);
      double angle = Mth::random(0.0, M_PI);
      double dx = Mth::random(-0.2, 0.2);
      double dy = Mth::random(-0.2, 0.2);

      OBB_2D obb(dx + centerX, dy + centerY, width, height, angle);

      obbs.push_back(obb);
    }
  }
  std::cout << "nb OBBs: " << obbs.size() << std::endl;
  size_t objectId = obbs.size()/2 - 1;

  //std::cout << "test 99-100: " << obbs[99].intersects(obbs[100]) << "\n";
  //std::cout << "test 100-101: " << obbs[100].intersects(obbs[101]) << "\n";

  //return 0;

  // O(N^2) complexity
  std::vector<std::pair<size_t, size_t>> VerletListSlow;
  std::vector<OBB_2D> intersOBBs1;
  auto start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < obbs.size(); i++) {
    for (size_t j=i+1; j < obbs.size(); j++) {
      //if (i == j)
      //  continue;
      if (obbs[i].intersects(obbs[j])) {
        VerletListSlow.push_back(std::make_pair(i, j));
        if (i == objectId)
          intersOBBs1.push_back(obbs[j]);
      }
    }
  }
  auto end = std::chrono::steady_clock::now();
  std::cout << "O(N^2) duration: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << " micro-seconds" << std::endl;

  // Generate SVG file with colored OBBs
  generateSVG("testOBB2D.svg", obbs[objectId], obbs, intersOBBs1, -18, -18, 18, 18);

  return 0;
}

#endif

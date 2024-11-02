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

#ifndef SVGTOOLS_HPP
#define SVGTOOLS_HPP
// https://www.alsacreations.com/tuto/lire/1421-svg-initiation-syntaxe-outils.html

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "vec2.hpp"

struct viewZone {
  double x, y, w, h;
  double scalex, scaley, x0, y0;

  viewZone(double x_, double y_, double w_, double h_) { setZone(x_, y_, w_, h_); }

  /// Set the zone of the view.
  ///
  /// The view is a rectangle defined by its lower left corner (x, y) and its
  /// width and height. The origin of the SVG coordinate system is at the lower
  /// left corner of the SVG file. The x-axis points to the right and the
  /// y-axis points upwards.
  void setZone(double x_, double y_, double w_, double h_) {
    x = x_;
    y = y_;
    w = w_;
    h = h_;
  }

  /// Adjust the range of the view.
  ///
  /// The range of the view is the bounding box of the data that is displayed.
  /// The x-axis and y-axis have different scales, so that the range is a
  /// rectangle on the SVG canvas. The origin of the SVG coordinate system is
  /// at the lower left corner of the SVG file. The x-axis points to the right
  /// and the y-axis points upwards.
  void adjustRange(double xmin, double xmax, double ymin, double ymax) {
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    scalex = w / dx;
    scaley = -h / dy;
    x0 = x - xmin * scalex;
    y0 = y + h - ymin * scaley;
  }
};

class SVGfile {
private:
  std::ostream &os;

public:
  SVGfile(std::ostream &outstream = std::cout) : os(outstream) {}

  /// @brief Initialize the SVG file with the given dimensions.
  ///
  /// This function writes the opening tags of an SVG file to the output stream,
  /// setting up the XML and SVG namespaces, version, and dimensions of the SVG canvas.
  ///
  /// @param width The width of the SVG canvas.
  /// @param height The height of the SVG canvas.
  void begin(int width, int height) {
    os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
       << "xmlns:xlink=\"http://www.w3.org/1999/xlink\""
       << " version=\"1.1\" width=\"" << width << "\" height=\"" << height << "\">\n";
  }

  /// @brief Write the closing tag of the SVG file.
  ///
  /// This function writes the "</svg>" tag to the output stream to close the SVG file.
  void end() { os << "</svg>\n"; }

  /// @brief Write a string to the SVG file.
  ///
  /// This function writes the given string followed by a newline character to the output stream.
  ///
  /// @param s The string to be written to the SVG file.
  void put(const char *s) { os << s << '\n'; }

  /// @brief Write a line to the SVG file.
  ///
  /// This function writes a line from point (x1, y1) to point (x2, y2) with the given style string to the SVG file.
  ///
  /// @param x1 The x-coordinate of the starting point of the line.
  /// @param y1 The y-coordinate of the starting point of the line.
  /// @param x2 The x-coordinate of the ending point of the line.
  /// @param y2 The y-coordinate of the ending point of the line.
  /// @param style The style string for the line, e.g. "stroke: black; stroke-width: 1.0".
  void line(double x1, double y1, double x2, double y2, const char *style) {
    os << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"" << style
       << "\" />\n";
  }

  /// @brief Write a line to the SVG file using the given view zone.
  ///
  /// This function writes a line from point (x1, y1) to point (x2, y2) with the given style string to the SVG file.
  /// The coordinates are transformed according to the given view zone.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param x1 The x-coordinate of the starting point of the line.
  /// @param y1 The y-coordinate of the starting point of the line.
  /// @param x2 The x-coordinate of the ending point of the line.
  /// @param y2 The y-coordinate of the ending point of the line.
  /// @param style The style string for the line, e.g. "stroke: black; stroke-width: 1.0".
  void line(viewZone &vz, double x1, double y1, double x2, double y2, const char *style) {
    line(x1 * vz.scalex + vz.x0, y1 * vz.scaley + vz.y0, x2 * vz.scalex + vz.x0, y2 * vz.scaley + vz.y0, style);
  }

  /// @brief Write a rectangle to the SVG file.
  ///
  /// This function writes a rectangle with upper left corner (x, y) and size (w, h) with the given style string
  /// to the SVG file.
  ///
  /// @param x The x-coordinate of the upper left corner of the rectangle.
  /// @param y The y-coordinate of the upper left corner of the rectangle.
  /// @param w The width of the rectangle.
  /// @param h The height of the rectangle.
  /// @param style The style string for the rectangle, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void rect(double x, double y, double w, double h, const char *style) {
    os << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << h << "\" style=\"" << style
       << "\" />\n";
  }

  /// @brief Write a rectangle to the SVG file using the given view zone.
  ///
  /// This function writes a rectangle with upper left corner (x, y) and size (w, h) with the given style string
  /// to the SVG file. The coordinates are transformed according to the given view zone.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param x The x-coordinate of the upper left corner of the rectangle.
  /// @param y The y-coordinate of the upper left corner of the rectangle.
  /// @param w The width of the rectangle.
  /// @param h The height of the rectangle.
  /// @param style The style string for the rectangle, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void rect(viewZone &vz, double x, double y, double w, double h, const char *style) {
    os << "<rect x=\"" << x * vz.scalex + vz.x0 << "\" y=\"" << y * vz.scaley + vz.y0 << "\" width=\""
       << fabs(w * vz.scalex) << "\" height=\"" << fabs(h * vz.scaley) << "\" style=\"" << style << "\" />\n";
  }

  /// @brief Write a polygon to the SVG file.
  ///
  /// This function writes a polygon with the given vertices and style string to the SVG file.
  ///
  /// @param vert The vertices of the polygon.
  /// @param style The style string for the polygon, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void polygon(std::vector<vec2r> &vert, const char *style) {
    os << "<polygon points=\"";
    for (size_t i = 0; i < vert.size(); i++) {
      os << vert[i].x << "," << vert[i].y << " ";
    }
    os << "\" style=\"" << style << "\" />\n";
  }

  /// @brief Write a polygon to the SVG file using the given view zone.
  ///
  /// This function writes a polygon with the given vertices and style string to the SVG file.
  /// The coordinates are transformed according to the given view zone.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param vert The vertices of the polygon.
  /// @param style The style string for the polygon, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void polygon(viewZone &vz, std::vector<vec2r> &vert, const char *style) {
    os << "<polygon points=\"";
    for (size_t i = 0; i < vert.size(); i++) {
      os << vert[i].x * vz.scalex + vz.x0 << "," << vert[i].y * vz.scaley + vz.y0 << " ";
    }
    os << "\" style=\"" << style << "\" />\n";
  }

  /// @brief Write a circle to the SVG file.
  ///
  /// This function writes a circle with center (cx, cy) and radius r to the SVG file
  /// using the specified style string.
  ///
  /// @param cx The x-coordinate of the center of the circle.
  /// @param cy The y-coordinate of the center of the circle.
  /// @param r The radius of the circle.
  /// @param style The style string for the circle, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void circle(double cx, double cy, double r, const char *style) {
    os << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r << "\" style=\"" << style << "\" />\n";
  }

  /// @brief Write a circle to the SVG file using the given view zone.
  ///
  /// This function writes a circle with center (cx, cy) and radius r to the SVG file
  /// using the specified style string. The coordinates are transformed according to the given view zone.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param cx The x-coordinate of the center of the circle.
  /// @param cy The y-coordinate of the center of the circle.
  /// @param r The radius of the circle.
  /// @param style The style string for the circle, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void circle(viewZone &vz, double cx, double cy, double r, const char *style) { ellipse(vz, cx, cy, r, r, style); }

  /// @brief Write an ellipse to the SVG file.
  ///
  /// This function writes an ellipse with center (cx, cy) and radii (rx, ry) to the SVG file
  /// using the specified style string.
  ///
  /// @param cx The x-coordinate of the center of the ellipse.
  /// @param cy The y-coordinate of the center of the ellipse.
  /// @param rx The x-radius of the ellipse.
  /// @param ry The y-radius of the ellipse.
  /// @param style The style string for the ellipse, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void ellipse(double cx, double cy, double rx, double ry, const char *style) {
    os << "<ellipse cx=\"" << cx << "\" cy=\"" << cy << "\" rx=\"" << rx << "\" ry=\"" << ry << "\" style=\"" << style
       << "\" />\n";
  }

  /// @brief Write an ellipse to the SVG file using the given view zone.
  ///
  /// This function writes an ellipse with center (cx, cy) and radii (rx, ry) to the SVG file
  /// using the specified style string. The coordinates are transformed according to the given view zone.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param cx The x-coordinate of the center of the ellipse.
  /// @param cy The y-coordinate of the center of the ellipse.
  /// @param rx The x-radius of the ellipse.
  /// @param ry The y-radius of the ellipse.
  /// @param style The style string for the ellipse, e.g. "stroke: black; stroke-width: 1.0; fill: none".
  void ellipse(viewZone &vz, double cx, double cy, double rx, double ry, const char *style) {
    os << "<ellipse cx=\"" << cx * vz.scalex + vz.x0 << "\" cy=\"" << cy * vz.scaley + vz.y0 << "\" rx=\""
       << fabs(rx * vz.scalex) << "\" ry=\"" << fabs(ry * vz.scaley) << "\" style=\"" << style << "\" />\n";
  }

  /// @brief Write text to the SVG file.
  ///
  /// This function writes a text element at the specified (x, y) position
  /// with the given style and content to the SVG file.
  ///
  /// @param x The x-coordinate of the text position.
  /// @param y The y-coordinate of the text position.
  /// @param style The style string for the text, e.g. "font-size: 12px; fill: black".
  /// @param content The textual content to be displayed.
  void text(double x, double y, const char *style, const char *content) {
    os << "<text x=\"" << x << "\" y=\"" << y << "\" style=\"" << style << "\">\n";
    os << content << '\n';
    os << "</text>\n";
  }

  /// @brief Write text to the SVG file using the given view zone.
  ///
  /// This function writes a text element at the transformed (x, y) position
  /// with the given style and content to the SVG file.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param x The x-coordinate of the text position.
  /// @param y The y-coordinate of the text position.
  /// @param style The style string for the text, e.g. "font-size: 12px; fill: black".
  /// @param content The textual content to be displayed.
  void text(viewZone &vz, double x, double y, const char *style, const char *content) {
    os << "<text x=\"" << x * vz.scalex + vz.x0 << "\" y=\"" << y * vz.scaley + vz.y0 << "\" style=\"" << style
       << "\">\n";
    os << content << '\n';
    os << "</text>\n";
  }

  /// @brief Begin a new group element in the SVG file.
  ///
  /// This function writes a group element with the specified id to the SVG file,
  /// allowing subsequent elements to be grouped together.
  ///
  /// @param id The identifier for the SVG group element.
  void group(const char *id) { os << "<g id=\"" << id << "\">\n"; }

  /// @brief End the current group element in the SVG file.
  ///
  /// This function writes the closing tag for the most recently opened group
  /// element in the SVG file, effectively marking the end of that group.
  void endGroup() { os << "</g>\n"; }

  /// @brief Write a use element to the SVG file.
  ///
  /// This function writes a use element with the specified id, x and y coordinates
  /// to the SVG file. The use element is used to reference another element in the
  /// current SVG file and display it at the specified position.
  ///
  /// @param id The identifier of the element to reference.
  /// @param x The x-coordinate of the element to display.
  /// @param y The y-coordinate of the element to display.
  void use(const char *id, double x, double y) {
    os << "<use xlink:href=\"#" << id << "\" x=\"" << x << "\" y=\"" << y << "\"/>\n";
  }

  /// @brief Write a curve to the SVG file using the given view zone.
  ///
  /// This function writes a path element to the SVG file, representing a curve
  /// that passes through the given x and y coordinates. The coordinates are
  /// transformed according to the specified view zone.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param x The x-coordinates of the points through which the curve passes.
  /// @param y The y-coordinates of the points through which the curve passes.
  /// @param style The style string for the path, e.g. "stroke: black; fill: none".
  /// @note The x and y vectors must have the same size; otherwise, an error message
  /// will be printed and the function will return without writing anything.
  void curve(viewZone &vz, std::vector<double> &x, std::vector<double> &y, const char *style) {
    if (x.size() != y.size()) {
      std::cerr << "Error: x and y vectors must have the same size.\n";
      return;
    }

    os << "<path d=\"";
    os << "M " << vz.scalex * x[0] + vz.x0 << "," << vz.scaley * y[0] + vz.y0;
    for (size_t i = 1; i < x.size(); i++) {
      os << " L " << vz.scalex * x[i] + vz.x0 << "," << vz.scaley * y[i] + vz.y0;
    }
    os << "\" style=\"" << style << "\" />\n";
  }

  /// @brief Write a scatter plot to the SVG file using the given view zone.
  ///
  /// This function writes a collection of circle elements to the SVG file,
  /// representing a scatter plot. The x and y coordinates of the points are
  /// transformed according to the specified view zone. All circles are drawn
  /// with the same radius and style.
  ///
  /// @param vz The view zone to be used for the transformation.
  /// @param x The x-coordinates of the points to be plotted.
  /// @param y The y-coordinates of the points to be plotted.
  /// @param radius The radius of each circle.
  /// @param style The style string for the circles, e.g. "stroke: black; fill: none".
  /// @note The x and y vectors must have the same size; otherwise, an error message
  /// will be printed and the function will return without writing anything.
  void scatter(viewZone &vz, std::vector<double> &x, std::vector<double> &y, double radius, const char *style) {
    if (x.size() != y.size()) {
      std::cerr << "Error: x and y vectors must have the same size.\n";
      return;
    }

    for (size_t i = 0; i < x.size(); i++) {
      double cx = vz.scalex * x[i] + vz.x0;
      double cy = vz.scaley * y[i] + vz.y0;
      os << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << radius << "\" style=\"" << style << "\" />\n";
    }
  }
};

#endif /* end of include guard: SVGTOOLS_HPP */

#if 0
int main(int argc, char const *argv[]) {
  std::ofstream ofs("test.svg");
  SVGfile svg(ofs);

  svg.begin(200, 200);

  svg.put("<defs>");
  svg.group("box");
  svg.line(10, 10, 90, 90, "stroke:red");
  svg.rect(30, 30, 40, 40, "stroke:blue;fill:none");
  svg.circle(50, 50, 10, "stroke:green;fill:pink");
  svg.endGroup();
  svg.put("</defs>");
  svg.use("box", 0, 0);
  svg.use("box", 100, 0);
  svg.use("box", 0, 100);
  svg.use("box", 100, 100);

  std::vector<double> xvec, yvec;
  for (double x = -10; x <= 10; x += 0.2) {
    xvec.push_back(x);
    yvec.push_back(cos(x));
  }
  viewZone vz(0, 0, 100, 100);
  vz.adjustRange(-10, 10, -1.1, 1.1);
  svg.curve(vz, xvec, yvec, "stroke:blue;fill:none");
  svg.scatter(vz, xvec, yvec,2, "stroke:blue;fill:none");
  svg.line(vz, 0, 0.5, 10, 0.5, "stroke:red;fill:none");

  svg.text(100, 100, "stroke:none;fill:black", "coucou");

  svg.end();

  // printf("%x", 0);

  return 0;
}
#endif

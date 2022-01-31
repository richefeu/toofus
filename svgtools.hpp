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

struct viewZone {
  double x, y, w, h;
  double scalex, scaley, x0, y0;

  viewZone(double x_, double y_, double w_, double h_) { setZone(x_, y_, w_, h_); }

  void setZone(double x_, double y_, double w_, double h_) {
    x = x_;
    y = y_;
    w = w_;
    h = h_;
  }

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

  void begin(int width, int height) {
    os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
       << "xmlns:xlink=\"http://www.w3.org/1999/xlink\""
       << " version=\"1.1\" width=\"" << width << "\" height=\"" << height << "\">\n";
  }

  void end() { os << "</svg>\n"; }

  void put(const char *s) { os << s << '\n'; }

  void line(double x1, double y1, double x2, double y2, const char *style) {
    os << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"" << style
       << "\" />\n";
  }

  void line(viewZone &vz, double x1, double y1, double x2, double y2, const char *style) {
    line(x1 * vz.scalex + vz.x0, y1 * vz.scaley + vz.y0, x2 * vz.scalex + vz.x0, y2 * vz.scaley + vz.y0, style);
  }

  void rect(double x, double y, double w, double h, const char *style) {
    os << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << h << "\" style=\"" << style
       << "\" />\n";
  }

  void circle(double cx, double cy, double r, const char *style) {
    os << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r << "\" style=\"" << style << "\" />\n";
  }
  
  void circle(viewZone &vz, double cx, double cy, double r, const char *style) {
    ellipse(vz, cx, cy, r, r, style);
  }

  void ellipse(double cx, double cy, double rx, double ry, const char *style) {
    os << "<ellipse cx=\"" << cx << "\" cy=\"" << cy << "\" rx=\"" << rx << "\" ry=\"" << ry << "\" style=\"" << style
       << "\" />\n";
  }

  void ellipse(viewZone &vz, double cx, double cy, double rx, double ry, const char *style) {
    os << "<ellipse cx=\"" << cx * vz.scalex + vz.x0 << "\" cy=\"" << cy * vz.scaley + vz.y0 << "\" rx=\""
       << fabs(rx * vz.scalex) << "\" ry=\"" << fabs(ry * vz.scaley) << "\" style=\"" << style << "\" />\n";
  }

  void text(double x, double y, const char *style, const char *content) {
    os << "<text x=\"" << x << "\" y=\"" << y << "\" style=\"" << style << "\">\n";
    os << content << '\n';
    os << "</text>\n";
  }

  void text(viewZone &vz, double x, double y, const char *style, const char *content) {
    os << "<text x=\"" << x * vz.scalex + vz.x0 << "\" y=\"" << y * vz.scaley + vz.y0 << "\" style=\"" << style
       << "\">\n";
    os << content << '\n';
    os << "</text>\n";
  }

  void group(const char *id) { os << "<g id=\"" << id << "\">\n"; }

  void endGroup() { os << "</g>\n"; }

  void use(const char *id, double x, double y) {
    os << "<use xlink:href=\"#" << id << "\" x=\"" << x << "\" y=\"" << y << "\"/>\n";
  }

  void curve(viewZone &vz, std::vector<double> &x, std::vector<double> &y, const char *style) {
    os << "<path d=\"";
    os << "M " << vz.scalex * x[0] + vz.x0 << "," << vz.scaley * y[0] + vz.y0;
    for (size_t i = 1; i < x.size(); i++) {
      os << " L " << vz.scalex * x[i] + vz.x0 << "," << vz.scaley * y[i] + vz.y0;
    }
    os << "\" style=\"" << style << "\" />\n";
  }
};

#endif /* end of include guard: SVGTOOLS_HPP */

#if 0
int main(int argc, char const *argv[]) {
  std::ofstream ofs("test.svg");
  SVGfile svg(ofs);
  // SVGfile svg;

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
  svg.line(vz, 0, 0.5, 10, 0.5, "stroke:red;fill:none");

  svg.text(100, 100, "stroke:none;fill:black", "coucou");

  svg.end();

  // printf("%x", 0);

  return 0;
}
#endif

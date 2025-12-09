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

#ifndef GLTOOLS_HPP
#define GLTOOLS_HPP

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#include <cmath>
#include <cstdarg> // for va_start and va_end
#include <cstring> // for strcpy

#include <cstdlib>
#include <functional>

#include "ColorTable.hpp"
#include "OBB.hpp"

class facetSphere {
protected:
  static void normalize(GLfloat *a) {
    GLfloat d = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= d;
    a[1] /= d;
    a[2] /= d;
  }

  static void drawtri(GLfloat *a, GLfloat *b, GLfloat *c, int div, float r) {
    if (div <= 0) {
      glNormal3fv(a);
      glVertex3f(-a[0] * r, -a[1] * r, -a[2] * r);
      glNormal3fv(b);
      glVertex3f(-b[0] * r, -b[1] * r, -b[2] * r);
      glNormal3fv(c);
      glVertex3f(-c[0] * r, -c[1] * r, -c[2] * r);
    } else {
      GLfloat ab[3], ac[3], bc[3];
      for (int i = 0; i < 3; i++) {
        ab[i] = 0.5f * (a[i] + b[i]);
        ac[i] = 0.5f * (a[i] + c[i]);
        bc[i] = 0.5f * (b[i] + c[i]);
      }
      normalize(ab);
      normalize(ac);
      normalize(bc);
      drawtri(a, ab, ac, div - 1, r);
      drawtri(b, bc, ab, div - 1, r);
      drawtri(c, ac, bc, div - 1, r);
      drawtri(ab, bc, ac, div - 1, r);
    }
  }

public:
  static GLfloat vdata[12][3];
  static GLuint tindices[20][3];

  static void draw(int ndiv, float radius) {
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < 20; ++i) {
      drawtri(vdata[tindices[i][0]], vdata[tindices[i][1]], vdata[tindices[i][2]], ndiv, radius);
    }
    glEnd();
  }
};

#define X 0.525731112119133606f
#define Z 0.850650808352039932f
GLfloat facetSphere::vdata[12][3] = {{-X, 0.0f, Z}, {X, 0.0f, Z},  {-X, 0.0f, -Z}, {X, 0.0f, -Z},
                                     {0.0f, Z, X},  {0.0f, Z, -X}, {0.0f, -Z, X},  {0.0f, -Z, -X},
                                     {Z, X, 0.0f},  {-Z, X, 0.0f}, {Z, -X, 0.0f},  {-Z, -X, 0.0f}};
#undef X
#undef Z

GLuint facetSphere::tindices[20][3] = {{0, 4, 1},  {0, 9, 4},  {9, 5, 4},  {4, 5, 8},  {4, 8, 1},
                                       {8, 10, 1}, {8, 3, 10}, {5, 3, 8},  {5, 2, 3},  {2, 7, 3},
                                       {7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6}, {0, 1, 6},
                                       {6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5},  {7, 2, 11}};

class glShape {
public:
  static void sphere(float radius, int ndiv) {
    facetSphere::draw(ndiv, radius);
  }

  static void arrow(const vec3r &orig, const vec3r &arrow, double arrowSize = -1.0, double arrowAngle = 0.7) {
    vec3r dest = orig + arrow;

    glLineWidth(2.0f);
    glBegin(GL_LINES);
    glVertex3d(orig.x, orig.y, orig.z);
    glVertex3d(dest.x, dest.y, dest.z);
    glEnd();

    vec3r v    = arrow;
    double len = v.normalize();
    if (arrowSize <= 0.0) { arrowSize = 0.04 * len; }
    vec3r vmz(v.x, v.y, v.z - 1.0); // v - z

    vec3r a;
    if (norm2(vmz) > 0.1) {
      a.set(v.y, -v.x, 0.0);
    } else {
      a.set(-v.z, 0.0, v.x);
    }
    a.normalize();
    vec3r b = cross(v, a);

    vec3r head = dest - arrowSize * v;
    vec3r c;
    double r = arrowSize * tan(0.5 * arrowAngle);
    glBegin(GL_TRIANGLE_FAN);
    glVertex3d(dest.x, dest.y, dest.z);
    for (double angle = 0.0; angle <= 2.0 * M_PI; angle += 0.2 * M_PI) {
      c = cos(angle) * a + sin(angle) * b;
      glNormal3d(c.x, c.y, c.z); // Pas tout à fait juste (!)
      c = head + r * c;
      glVertex3d(c.x, c.y, c.z);
    }
    glEnd();
  }

  static void tube(vec3r &orig, vec3r &arrow, double diam) {
    vec3r dest = orig + arrow;
    vec3r v    = arrow;
    v.normalize();
    vec3r vmz(v.x, v.y, v.z - 1.0); // v - z

    vec3r a;
    if (norm2(vmz) > 0.1) {
      a.set(v.y, -v.x, 0.0);
    } else {
      a.set(-v.z, 0.0, v.x);
    }
    if (a.isnull()) {
      a.set(1.0, 1.0, 0.0); // for the rare cases v = (0, 0, -1)
    }
    a.normalize();
    vec3r b = cross(v, a);

    vec3r c1, c2, n;
    double r = 0.5 * diam;
    glBegin(GL_TRIANGLE_STRIP);
    for (double angle = 0.0; angle <= 2.0 * M_PI; angle += 0.05 * M_PI) {
      n = cos(angle) * a + sin(angle) * b;
      glNormal3d(n.x, n.y, n.z);
      n *= r;
      c1 = orig + n;
      c2 = dest + n;
      glVertex3d(c1.x, c1.y, c1.z);
      glVertex3d(c2.x, c2.y, c2.z);
    }
    glEnd();
  }

  // Attention : jamais testé !!!!!
  static void triangle(vec3r &A, vec3r &B, vec3r &C, vec3r &NA, vec3r &NB, vec3r &NC) {
    glNormal3d(NA.x, NA.y, NA.z);
    glVertex3d(A.x, A.y, A.z);
    glNormal3d(NB.x, NB.y, NB.z);
    glVertex3d(B.x, B.y, B.z);
    glNormal3d(NC.x, NC.y, NC.z);
    glVertex3d(C.x, C.y, C.z);
  }

  static void frame(const vec3r &pos, double lx = 1.0, double ly = 1.0, double lz = 1.0) {
    glDisable(GL_LIGHTING);

    double arrowSize = lx;
    if (ly > arrowSize) { arrowSize = ly; }
    if (lz > arrowSize) { arrowSize = lz; }
    arrowSize *= 0.02;

    glColor4ub(255, 0, 0, 255);
    glShape::arrow(pos, lx * vec3r::unit_x(), arrowSize);

    glColor4ub(0, 255, 0, 255);
    glShape::arrow(pos, ly * vec3r::unit_y(), arrowSize);

    glColor4ub(0, 0, 255, 255);
    glShape::arrow(pos, lz * vec3r::unit_z(), arrowSize);
  }

  static void obb(OBB &obb) {
    glDisable(GL_LIGHTING);

    glLineWidth(1.0f);

    vec3r orig = obb.center;
    vec3r e0   = obb.extent[0] * obb.e[0];
    vec3r e1   = obb.extent[1] * obb.e[1];
    vec3r e2   = obb.extent[2] * obb.e[2];

    vec3r p0 = orig - e0 + e1 + e2;
    vec3r p1 = orig - e0 - e1 + e2;
    vec3r p2 = orig - e0 - e1 - e2;
    vec3r p3 = orig - e0 + e1 - e2;

    vec3r p4 = orig + e0 + e1 + e2;
    vec3r p5 = orig + e0 - e1 + e2;
    vec3r p6 = orig + e0 - e1 - e2;
    vec3r p7 = orig + e0 + e1 - e2;

    glBegin(GL_LINE_LOOP);
    glVertex3d(p0.x, p0.y, p0.z);
    glVertex3d(p1.x, p1.y, p1.z);
    glVertex3d(p2.x, p2.y, p2.z);
    glVertex3d(p3.x, p3.y, p3.z);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex3d(p4.x, p4.y, p4.z);
    glVertex3d(p5.x, p5.y, p5.z);
    glVertex3d(p6.x, p6.y, p6.z);
    glVertex3d(p7.x, p7.y, p7.z);
    glEnd();

    glBegin(GL_LINES);
    glVertex3d(p0.x, p0.y, p0.z);
    glVertex3d(p4.x, p4.y, p4.z);

    glVertex3d(p1.x, p1.y, p1.z);
    glVertex3d(p5.x, p5.y, p5.z);

    glVertex3d(p2.x, p2.y, p2.z);
    glVertex3d(p6.x, p6.y, p6.z);

    glVertex3d(p3.x, p3.y, p3.z);
    glVertex3d(p7.x, p7.y, p7.z);
    glEnd();
  }
};

// Usage:
// switch2D::go(w, h);
// Draw things with the coordinate (0, 0) sets at the bottom left of the window
// switch2D::back(); // All attributes and matrix are set back
class switch2D {
public:
  static void go(int w, int h) {
    // (x,y) is from the bottom left of the window
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, w, 0, h, -1.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glPushAttrib(GL_DEPTH_TEST);
    glDisable(GL_DEPTH_TEST);
    glPushAttrib(GL_LIGHTING);
    glDisable(GL_LIGHTING);
  }

  static void back() {
    glPopAttrib();
    glPopAttrib();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
  }
};

class glText {
public:
  static void makeRasterFont() {
    GLuint i;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    fontOffset = glGenLists(128);
    for (i = 32; i < 127; i++) {
      glNewList(i + fontOffset, GL_COMPILE);
      glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, rasters[i - 32]);
      glEndList();
    }
  }

  static GLubyte rasters[95][13];
  static GLuint fontOffset;

  // DO NOT FORGET TO INIT BEFORE BEING ABLE TO PRINT
  static void init() {
    glShadeModel(GL_FLAT);
    makeRasterFont();
  }

  static void print(int x, int y, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char buffer[128];
    vsnprintf(buffer, 127, fmt, args);
    va_end(args);

    glRasterPos2i(x, y);
    glPushAttrib(GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists((GLsizei)strlen(buffer), GL_UNSIGNED_BYTE, (GLubyte *)buffer);
    glPopAttrib();
  }

  static void print(GLfloat x, GLfloat y, GLfloat z, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char buffer[128];
    vsnprintf(buffer, 127, fmt, args);
    va_end(args);

    glRasterPos3f(x, y, z);
    glPushAttrib(GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists((GLsizei)strlen(buffer), GL_UNSIGNED_BYTE, (GLubyte *)buffer);
    glPopAttrib();
  }
};

GLubyte glText::rasters[95][13] = {{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18},
                                   {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x36, 0x36, 0x36},
                                   {0x00, 0x00, 0x00, 0x66, 0x66, 0xff, 0x66, 0x66, 0xff, 0x66, 0x66, 0x00, 0x00},
                                   {0x00, 0x00, 0x18, 0x7e, 0xff, 0x1b, 0x1f, 0x7e, 0xf8, 0xd8, 0xff, 0x7e, 0x18},
                                   {0x00, 0x00, 0x0e, 0x1b, 0xdb, 0x6e, 0x30, 0x18, 0x0c, 0x76, 0xdb, 0xd8, 0x70},
                                   {0x00, 0x00, 0x7f, 0xc6, 0xcf, 0xd8, 0x70, 0x70, 0xd8, 0xcc, 0xcc, 0x6c, 0x38},
                                   {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x1c, 0x0c, 0x0e},
                                   {0x00, 0x00, 0x0c, 0x18, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c},
                                   {0x00, 0x00, 0x30, 0x18, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x18, 0x30},
                                   {0x00, 0x00, 0x00, 0x00, 0x99, 0x5a, 0x3c, 0xff, 0x3c, 0x5a, 0x99, 0x00, 0x00},
                                   {0x00, 0x00, 0x00, 0x18, 0x18, 0x18, 0xff, 0xff, 0x18, 0x18, 0x18, 0x00, 0x00},
                                   {0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x60, 0x60, 0x30, 0x30, 0x18, 0x18, 0x0c, 0x0c, 0x06, 0x06, 0x03, 0x03},
                                   {0x00, 0x00, 0x3c, 0x66, 0xc3, 0xe3, 0xf3, 0xdb, 0xcf, 0xc7, 0xc3, 0x66, 0x3c},
                                   {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78, 0x38, 0x18},
                                   {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0xe7, 0x7e},
                                   {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0x07, 0x03, 0x03, 0xe7, 0x7e},
                                   {0x00, 0x00, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0xff, 0xcc, 0x6c, 0x3c, 0x1c, 0x0c},
                                   {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
                                   {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
                                   {0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x03, 0x03, 0xff},
                                   {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e},
                                   {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x03, 0x7f, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e},
                                   {0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x1c, 0x1c, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x06, 0x0c, 0x18, 0x30, 0x60, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06},
                                   {0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x06, 0x0c, 0x18, 0x30, 0x60},
                                   {0x00, 0x00, 0x18, 0x00, 0x00, 0x18, 0x18, 0x0c, 0x06, 0x03, 0xc3, 0xc3, 0x7e},
                                   {0x00, 0x00, 0x3f, 0x60, 0xcf, 0xdb, 0xd3, 0xdd, 0xc3, 0x7e, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18},
                                   {0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
                                   {0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
                                   {0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc},
                                   {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
                                   {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff},
                                   {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
                                   {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
                                   {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e},
                                   {0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06},
                                   {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3},
                                   {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
                                   {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3},
                                   {0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3},
                                   {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e},
                                   {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
                                   {0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c},
                                   {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
                                   {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e},
                                   {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff},
                                   {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
                                   {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
                                   {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
                                   {0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
                                   {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
                                   {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff},
                                   {0x00, 0x00, 0x3c, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x3c},
                                   {0x00, 0x03, 0x03, 0x06, 0x06, 0x0c, 0x0c, 0x18, 0x18, 0x30, 0x30, 0x60, 0x60},
                                   {0x00, 0x00, 0x3c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x3c},
                                   {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18},
                                   {0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x38, 0x30, 0x70},
                                   {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0x7f, 0x03, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
                                   {0x00, 0x00, 0x7e, 0xc3, 0xc0, 0xc0, 0xc0, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x03, 0x03, 0x03, 0x03},
                                   {0x00, 0x00, 0x7f, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x33, 0x1e},
                                   {0x7e, 0xc3, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0},
                                   {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x00, 0x18, 0x00},
                                   {0x38, 0x6c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x00, 0x00, 0x0c, 0x00},
                                   {0x00, 0x00, 0xc6, 0xcc, 0xf8, 0xf0, 0xd8, 0xcc, 0xc6, 0xc0, 0xc0, 0xc0, 0xc0},
                                   {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78},
                                   {0x00, 0x00, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xfe, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xfc, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x7c, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x7c, 0x00, 0x00, 0x00, 0x00},
                                   {0xc0, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00},
                                   {0x03, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe0, 0xfe, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xfe, 0x03, 0x03, 0x7e, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x1c, 0x36, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x00},
                                   {0x00, 0x00, 0x7e, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xdb, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18, 0x3c, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00},
                                   {0xc0, 0x60, 0x60, 0x30, 0x18, 0x3c, 0x66, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0xff, 0x60, 0x30, 0x18, 0x0c, 0x06, 0xff, 0x00, 0x00, 0x00, 0x00},
                                   {0x00, 0x00, 0x0f, 0x18, 0x18, 0x18, 0x38, 0xf0, 0x38, 0x18, 0x18, 0x18, 0x0f},
                                   {0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18},
                                   {0x00, 0x00, 0xf0, 0x18, 0x18, 0x18, 0x1c, 0x0f, 0x1c, 0x18, 0x18, 0x18, 0xf0},
                                   {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x8f, 0xf1, 0x60, 0x00, 0x00, 0x00}};

GLuint glText::fontOffset = 0;

class glTextZone {
#define NB_LINE_MAX 40

protected:
  char textzone[NB_LINE_MAX][128];
  int nbLine;
  int *width;
  int *height;

  void drawBackground() {
    glColor4f(1.0f, 1.0f, 1.0f, 0.6f);
    int h = nbLine * 16;
    glBegin(GL_QUADS);
    glVertex2i(0, 0);
    glVertex2i(*width, 0);
    glVertex2i(*width, h);
    glVertex2i(0, h);
    glEnd();
    
    glColor4f(0.0f, 0.0f, 0.0f, 0.6f);
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    glVertex2i(0, h);
    glVertex2i(*width, h);
    glEnd();
  }

public:
  // Ctors
  glTextZone(int *W, int *H) : nbLine(NB_LINE_MAX), width(W), height(H) {}
  glTextZone(int n, int *W, int *H) : nbLine(n), width(W), height(H) {}

  void set_nbLine(int nb) {
    if (nb > 0 && nb < NB_LINE_MAX) { nbLine = nb; }
  }

  void set_min_nbLine(int min) {
    if (nbLine < min && min > 0 && min < NB_LINE_MAX) { nbLine = min; }
  }

  void increase_nbLine() {
    nbLine++;
    if (nbLine >= NB_LINE_MAX) { nbLine = NB_LINE_MAX - 1; }
  }

  void decrease_nbLine() {
    nbLine--;
    if (nbLine < 1) { nbLine = 1; }
  }

  void reset() {
    for (size_t i = 0; i < NB_LINE_MAX; i++) { textzone[i][0] = '\0'; }
  }

  void draw() {
    switch2D::go(*width, *height);
    drawBackground();

    glColor3i(0, 0, 0);
    for (int i = 0; i < nbLine; ++i) { glText::print(4, 4 + i * 16, textzone[i]); }

    switch2D::back();
  }

  void addLine(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char buffer[128];
    vsnprintf(buffer, 128, fmt, args);
    for (int i = NB_LINE_MAX - 1; i > 0; i--) { strncpy(textzone[i], textzone[i - 1], 128); }
    snprintf((char *)textzone[0], 128, "%s", buffer);
    va_end(args);
  }

#undef NB_LINE_MAX
};

class glColorBar {
protected:
  // position is from top left
  int xpos;
  int ypos;
  int wbox;
  int hbox;
  char title[128];
  std::vector<int> label_index_pos;
  std::vector<std::string> labels;

public:
  glColorBar() : xpos(10), ypos(24), wbox(25), hbox(200) {
    snprintf(title, 128, "notitle");
  }

  void addLabel(int i, std::string &lab, ColorTable &ct) {
    if (i >= 0 && i < ct.getSize()) {
      label_index_pos.push_back(i);
      labels.push_back(lab);
    }
  }

  void setPos(int x, int y) {
    xpos = x;
    ypos = y;
  }

  void setSize(int w, int h) {
    wbox = w;
    hbox = h;
  }

  void setTitle(const char *t) {
    strncpy(title, t, 128);
  }

  void show(int W, int H, ColorTable &ct) {
    switch2D::go(W, H);
    glLineWidth(1.0f);

    colorRGBA col;
    float value;
    float dval = (ct.getMax() - ct.getMin()) / (float)(ct.getSize() - 1);
    // std::cout << "dval = " << dval << std::endl;
    float dH     = (float)hbox / (float)(ct.getSize());
    float bottom = (float)(H - (ypos + hbox));

    // draw the color bar
    for (int i = 0; i < ct.getSize(); ++i) {
      value = ct.getMin() + (float)i * dval;
      // std::cout << "value = " << value << std::endl;
      ct.getRGB(value, &col);
      glColor3ub((GLubyte)col.r, (GLubyte)col.g, (GLubyte)col.b);
      glBegin(GL_QUADS);
      glVertex2i(xpos, (int)floor(bottom + (float)i * dH));
      glVertex2i((xpos + wbox), (int)floor(bottom + (float)i * dH));
      glVertex2i((xpos + wbox), (int)floor(bottom + (float)(i + 1) * dH));
      glVertex2i(xpos, (int)floor(bottom + (float)(i + 1) * dH));
      glEnd();
    }

    // surrounding box
    glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(xpos, H - ypos);
    glVertex2i(xpos + wbox, H - ypos);
    glVertex2i(xpos + wbox, H - (ypos + hbox));
    glVertex2i(xpos, H - (ypos + hbox));
    glEnd();

    if (labels.empty()) {
      char strValMin[64];
      char strValMax[64];
      snprintf(strValMin, 64, "%.1E", ct.getMin());
      snprintf(strValMax, 64, "%.1E", ct.getMax());
      glText::print(xpos + wbox + 4, (int)floor(bottom), strValMin);
      glText::print(xpos + wbox + 4, (int)floor(bottom) + hbox - 13, strValMax);

    } else {
      for (size_t i = 0; i < labels.size(); i++) {
        glText::print(xpos + wbox + 4, (int)floor(bottom + (float)(label_index_pos[i] + 0.5) * dH - 7),
                      labels[i].c_str());
      }
    }

    glText::print(xpos, H - ypos + 6, title);

    switch2D::back();
  }
};

class glTools {
public:
  // This funtion can be called before drawing anything.
  // It clears the screen with eventually a color gradient from bottom to top.
  // Default is light blue to white.
  static void clearBackground(bool grad, int bottomRed = 135, int bottomGreen = 206, int bottomBlue = 250,
                              int topRed = 255, int topGreen = 255, int topBlue = 255) {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (!grad) { return; }

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_QUADS);
    glColor3ub((GLubyte)bottomRed, (GLubyte)bottomGreen, (GLubyte)bottomBlue); // Bottom color
    glVertex2f(-1.0f, -1.0f);
    glVertex2f(1.0f, -1.0f);
    glColor3ub((GLubyte)topRed, (GLubyte)topGreen, (GLubyte)topBlue); // Top color
    glVertex2f(1.0f, 1.0f);
    glVertex2f(-1.0f, 1.0f);
    glEnd();

    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
  }
};

#endif /* end of include guard: GLTOOLS_HPP */

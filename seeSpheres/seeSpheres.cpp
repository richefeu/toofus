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

#include "seeSpheres.hpp"

void clear_background() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (!show_background) return;

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_QUADS);
  glColor3f(0.4f, 0.4f, 1.0f);  // Bottom color
  glVertex2f(-1.0f, -1.0f);
  glVertex2f(1.0f, -1.0f);
  glColor3f(1.0f, 1.0f, 1.0f);  // Top color
  glVertex2f(1.0f, 1.0f);
  glVertex2f(-1.0f, 1.0f);
  glEnd();

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
}

void keyboard(unsigned char Key, int x, int y) {
  switch (Key) {
    case 'b':
      show_background = 1 - show_background;
      break;

    case 'c': {
      rminRateConnection -= 0.1;
      rminRateConnection = std::max(0.0, rminRateConnection);
      buildPeriodicNeighbors();
    } break;

    case 'C': {
      rminRateConnection += 0.1;
      rminRateConnection = std::min(1.0, rminRateConnection);
      buildPeriodicNeighbors();
    } break;

    case 'q':
      exit(0);
      break;

    case 't': {
      transparency -= 0.1;
      transparency = std::max(0.0, transparency);
    } break;

    case 'T': {
      transparency += 0.1;
      transparency = std::min(1.0, transparency);
    } break;

    case '=': {
      fit_view();
      adjust_clipping_plans();
    } break;
  };

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
          mouse_mode = PAN;
        else
          mouse_mode = ROTATION;
        break;
      case GLUT_MIDDLE_BUTTON:
        mouse_mode = ZOOM;
        break;
    }
  }
}

void rotatePoint(double p_x, double p_y, double p_z, double center_x, double center_y, double center_z, double axis_x,
                 double axis_y, double axis_z, double& ret_x, double& ret_y, double& ret_z, double theta) {
  double c = cos(theta), s = sin(theta);
  double C = 1.0 - c;

  double tmp_x = p_x - center_x;
  double tmp_y = p_y - center_y;
  double tmp_z = p_z - center_z;

  ret_x = center_x + tmp_x * (axis_x * axis_x * C + c) + tmp_y * (axis_x * axis_y * C - axis_z * s) +
          tmp_z * (axis_x * axis_z * C + axis_y * s);
  ret_y = center_y + tmp_x * (axis_y * axis_x * C + axis_z * s) + tmp_y * (axis_y * axis_y * C + c) +
          tmp_z * (axis_y * axis_z * C - axis_x * s);
  ret_z = center_z + tmp_x * (axis_z * axis_x * C - axis_y * s) + tmp_y * (axis_z * axis_y * C + axis_x * s) +
          tmp_z * (axis_z * axis_z * C + c);
}

void motion(int x, int y) {
  if (mouse_mode == NOTHING) return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

    case ROTATION: {
      double center_eye_x = center_x - eye_x;
      double center_eye_y = center_y - eye_y;
      double center_eye_z = center_z - eye_z;
      double axis_x = up_y * center_eye_z - up_z * center_eye_y;
      double axis_y = up_z * center_eye_x - up_x * center_eye_z;
      double axis_z = up_x * center_eye_y - up_y * center_eye_x;
      double axis_len = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
      if (axis_len > 0.0) {
        axis_x /= axis_len;
        axis_y /= axis_len;
        axis_z /= axis_len;
      }
      double res_x = 0.0;
      double res_y = 0.0;
      double res_z = 0.0;
      rotatePoint(eye_x, eye_y, eye_z, center_x, center_y, center_z, up_x, up_y, up_z, res_x, res_y, res_z, -dx * M_PI);
      eye_x = res_x;
      eye_y = res_y;
      eye_z = res_z;
      rotatePoint(eye_x, eye_y, eye_z, center_x, center_y, center_z, axis_x, axis_y, axis_z, res_x, res_y, res_z,
                  dy * M_PI);
      eye_x = res_x;
      eye_y = res_y;
      eye_z = res_z;
      double center_up_x = center_x + up_x;
      double center_up_y = center_y + up_y;
      double center_up_z = center_z + up_z;
      rotatePoint(center_up_x, center_up_y, center_up_z, center_x, center_y, center_z, axis_x, axis_y, axis_z, res_x,
                  res_y, res_z, dy * M_PI);
      up_x = res_x - center_x;
      up_y = res_y - center_y;
      up_z = res_z - center_z;
      double up_len = sqrt(up_x * up_x + up_y * up_y + up_z * up_z);
      if (up_len > 0.0) {
        up_x /= up_len;
        up_y /= up_len;
        up_z /= up_len;
      }
    } break;

    case ZOOM: {
      eye_x = center_x + (eye_x - center_x) * (dy + 1.0);
      eye_y = center_y + (eye_y - center_y) * (dy + 1.0);
      eye_z = center_z + (eye_z - center_z) * (dy + 1.0);
    } break;

    case PAN: {
      double center_eye_x = center_x - eye_x;
      double center_eye_y = center_y - eye_y;
      double center_eye_z = center_z - eye_z;
      double center_eye_len =
          sqrt(center_eye_x * center_eye_x + center_eye_y * center_eye_y + center_eye_z * center_eye_z);

      double length = center_eye_len * tan(view_angle * M_PI / 360.0) * 2.0;
      double axis_x = up_y * center_eye_z - up_z * center_eye_y;
      double axis_y = up_z * center_eye_x - up_x * center_eye_z;
      double axis_z = up_x * center_eye_y - up_y * center_eye_x;
      double axis_len = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
      if (axis_len > 0.0) {
        axis_x /= axis_len;
        axis_y /= axis_len;
        axis_z /= axis_len;
      }
      center_x += (axis_x * dx * length * 0.8) + (up_x * dy * length);
      center_y += (axis_y * dx * length * 0.8) + (up_y * dy * length);
      center_z += (axis_z * dx * length * 0.8) + (up_z * dy * length);
    } break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  display();
}

void normalize(GLfloat* a) {
  GLfloat d = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  a[0] /= d;
  a[1] /= d;
  a[2] /= d;
}

void drawtri(GLfloat* a, GLfloat* b, GLfloat* c, int div, float r) {
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
      ab[i] = (a[i] + b[i]) / 2;
      ac[i] = (a[i] + c[i]) / 2;
      bc[i] = (b[i] + c[i]) / 2;
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

void drawsphere(int ndiv, float radius) {
  glBegin(GL_TRIANGLES);
  for (int i = 0; i < 20; ++i)
    drawtri(vdata[tindices[i][0]], vdata[tindices[i][1]], vdata[tindices[i][2]], ndiv, radius);
  glEnd();
}

void display() {
  clear_background();
  adjust_clipping_plans();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(eye_x, eye_y, eye_z, center_x, center_y, center_z, up_x, up_y, up_z);

  glShadeModel(GL_SMOOTH);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  // Display things
  drawParticles();
  drawConnections();
  drawBoundingBox();

  glFlush();
  glutSwapBuffers();
}

void adjust_clipping_plans() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  wh_ratio = (float)width / (float)height;
  double eye_center_x = eye_x - center_x;
  double eye_center_y = eye_y - center_y;
  double eye_center_z = eye_z - center_z;
  double zf = sqrt(eye_center_x * eye_center_x + eye_center_y * eye_center_y + eye_center_z * eye_center_z);
  double lx = xmax - xmin;
  double ly = ymax - ymin;
  double lz = zmax - zmin;
  max_length = sqrt(lx * lx + ly * ly + lz * lz);
  double znear = zf - 0.5 * max_length;
  double close_dst = 0.1 * zf;
  if (znear < close_dst) znear = close_dst;
  double zfar = zf + 0.5 * max_length;
  gluPerspective(view_angle, wh_ratio, znear, zfar);
  glMatrixMode(GL_MODELVIEW);
}

void fit_view() {
  double dir_x = eye_x - center_x;
  double dir_y = eye_y - center_y;
  double dir_z = eye_z - center_z;
  double diag_x = xmax - xmin;
  double diag_y = ymax - ymin;
  double diag_z = zmax - zmin;
  double dir_len = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
  if (dir_len > 0.0) {
    dir_x /= dir_len;
    dir_y /= dir_len;
    dir_z /= dir_len;
  }
  center_x = xmin + 0.5 * diag_x;
  center_y = ymin + 0.5 * diag_y;
  center_z = zmin + 0.5 * diag_z;
  double diag_len = sqrt(diag_x * diag_x + diag_y * diag_y + diag_z * diag_z);
  GLfloat d = 0.5 * diag_len / (atan(view_angle * M_PI / 360.0));
  eye_x = center_x + d * dir_x;
  eye_y = center_y + d * dir_y;
  eye_z = center_z + d * dir_z;
}

void reshape(int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, width, height);

  adjust_clipping_plans();
  glutPostRedisplay();
}

void drawTube(double origx, double origy, double origz, double arrowx, double arrowy, double arrowz, double diam) {
  double destx = origx + arrowx;
  double desty = origy + arrowy;
  double destz = origz + arrowz;

  double vx = arrowx;
  double vy = arrowy;
  double vz = arrowz;
  double l = sqrt(vx * vx + vy * vy + vz * vz);
  vx /= l;
  vy /= l;
  vz /= l;

  double vmzx = vx;
  double vmzy = vy;
  double vmzz = vz - 1.0;

  double ax, ay, az;
  double norm2_vmz = vmzx * vmzx + vmzy * vmzy + vmzz * vmzz;
  if (norm2_vmz > 0.1) {
    ax = vy;
    ay = -vx;
    az = 0.0;
  } else {
    ax = -vz;
    ay = 0.0;
    az = vx;
  }

  l = sqrt(ax * ax + ay * ay + az * az);
  ax /= l;
  ay /= l;
  az /= l;

  double bx = vy * az - vz * ay;
  double by = vz * ax - vx * az;
  double bz = vx * ay - vy * ax;

  double c1x, c1y, c1z;
  double c2x, c2y, c2z;
  double nx, ny, nz;
  double r = 0.5 * diam;
  glBegin(GL_TRIANGLE_STRIP);
  for (double angle = 0.0; angle <= 2.0 * M_PI; angle += 0.2 * M_PI) {
    double c = cos(angle);
    double s = sin(angle);
    nx = c * ax + s * bx;
    ny = c * ay + s * by;
    nz = c * az + s * bz;
    glNormal3f(nx, ny, nz);
    nx *= r;
    ny *= r;
    nz *= r;

    c1x = origx + nx;
    c1y = origy + ny;
    c1z = origz + nz;

    c2x = destx + nx;
    c2y = desty + ny;
    c2z = destz + nz;
    glVertex3f(c1x, c1y, c1z);
    glVertex3f(c2x, c2y, c2z);
  }
  glEnd();
}

void drawParticles() {
  if (mouse_mode != NOTHING && spheres.size() > 2000) return;

  glColor4f(0.337254901960784, 0.505882352941176, 0.768627450980392, transparency);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  for (size_t i = 0; i < spheres.size(); ++i) {
    glPushMatrix();
    glTranslatef(spheres[i].x, spheres[i].y, spheres[i].z);
    drawsphere(2, spheres[i].r);
    glPopMatrix();
  }
}

void drawConnections() {
  if (mouse_mode != NOTHING && spheres.size() > 2000) return;

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glLineWidth(3.0f);
  glColor3f(1.0f, 0.5f, 0.0f);

  glBegin(GL_LINES);
  for (size_t k = 0; k < neighbors.size(); ++k) {
    int i = neighbors[k].i;
    int j = neighbors[k].j;
    double l = neighbors[k].dst;
    glVertex3f(spheres[i].x, spheres[i].y, spheres[i].z);
    glVertex3f(spheres[i].x + l * neighbors[k].nx, spheres[i].y + l * neighbors[k].ny,
               spheres[i].z + l * neighbors[k].nz);

    glVertex3f(spheres[j].x, spheres[j].y, spheres[j].z);
    glVertex3f(spheres[j].x - l * neighbors[k].nx, spheres[j].y - l * neighbors[k].ny,
               spheres[j].z - l * neighbors[k].nz);
  }
  glEnd();
}

void drawBoundingBox() {
  glDisable(GL_LIGHTING);

  glLineWidth(3.0f);
  glColor3f(1.0f, 0.5f, 0.84f);

  glBegin(GL_LINE_LOOP);
  glVertex3f(xmin, ymin, zmin);
  glVertex3f(xmax, ymin, zmin);
  glVertex3f(xmax, ymax, zmin);
  glVertex3f(xmin, ymax, zmin);
  glVertex3f(xmin, ymin, zmin);
  glEnd();

  glBegin(GL_LINE_LOOP);
  glVertex3f(xmin, ymin, zmax);
  glVertex3f(xmax, ymin, zmax);
  glVertex3f(xmax, ymax, zmax);
  glVertex3f(xmin, ymax, zmax);
  glVertex3f(xmin, ymin, zmax);
  glEnd();

  glBegin(GL_LINES);
  glVertex3f(xmin, ymin, zmin);
  glVertex3f(xmin, ymin, zmax);

  glVertex3f(xmax, ymin, zmin);
  glVertex3f(xmax, ymin, zmax);

  glVertex3f(xmax, ymax, zmin);
  glVertex3f(xmax, ymax, zmax);

  glVertex3f(xmin, ymax, zmin);
  glVertex3f(xmin, ymax, zmax);
  glEnd();
}

// Robust and portable function to test if a file exists
bool fileExists(const char* fileName) {
  std::fstream fin;
  fin.open(fileName, std::ios::in);
  if (fin.is_open()) {
    fin.close();
    return true;
  }
  fin.close();
  return false;
}

void bounding_box() {
  xmin = ymin = zmin = 1e20;
  xmax = ymax = zmax = -1e20;
  for (size_t i = 0; i < spheres.size(); ++i) {
    if (spheres[i].x - spheres[i].r < xmin) xmin = spheres[i].x - spheres[i].r;
    if (spheres[i].y - spheres[i].r < ymin) ymin = spheres[i].y - spheres[i].r;
    if (spheres[i].z - spheres[i].r < zmin) zmin = spheres[i].z - spheres[i].r;
    if (spheres[i].x + spheres[i].r > xmax) xmax = spheres[i].x + spheres[i].r;
    if (spheres[i].y + spheres[i].r > ymax) ymax = spheres[i].y + spheres[i].r;
    if (spheres[i].z + spheres[i].r > zmax) zmax = spheres[i].z + spheres[i].r;
  }
}

void readSpheres(const char* name) {
  if (!fileExists(name)) {
    std::cout << name << " does not exist" << std::endl;
    return;
  }

  std::cout << "Read " << name << std::endl;
  std::ifstream is(name);

  while (is.good()) {
    if (is.eof()) break;
    Sphere S;
    is >> S.x >> S.y >> S.z >> S.r;
    if (is.eof()) break;
    spheres.push_back(S);

    if (is.eof()) break;
  }
  if (spheres.empty()) return;
  std::cout << "Number of spheres: " << spheres.size() << std::endl;

  if (fileExists("box.txt")) {
    std::ifstream boxfile("box.txt");
    boxfile >> xmin >> ymin >> zmin;
    boxfile >> xmax >> ymax >> zmax;
  } else {
    bounding_box();
  }
  buildPeriodicNeighbors();
}

void buildPeriodicNeighbors() {

  double rmin = 1e20;
  for (size_t i = 0; i < spheres.size(); i++) {
    rmin = std::min(rmin, spheres[i].r);
  }
  rmin *= rminRateConnection;

  neighbors.clear();
  double sx = xmax - xmin;
  double sy = ymax - ymin;
  double sz = zmax - zmin;

  for (size_t i = 0; i < spheres.size(); i++) {
    for (size_t j = i + 1; j < spheres.size(); j++) {
      double lijx = (spheres[j].x - spheres[i].x) / sx;
      double lijy = (spheres[j].y - spheres[i].y) / sy;
      double lijz = (spheres[j].z - spheres[i].z) / sz;
      lijx -= floor(lijx + 0.5);
      lijy -= floor(lijy + 0.5);
      lijz -= floor(lijz + 0.5);
      double bx = sx * lijx;
      double by = sy * lijy;
      double bz = sz * lijz;
      double ccdist = sqrt(bx * bx + by * by + bz * bz);
      if (ccdist < spheres[i].r + spheres[j].r + rmin) {
        Neighbor N;
        N.i = i;
        N.j = j;
        N.nx = bx / ccdist;
        N.ny = by / ccdist;
        N.nz = bz / ccdist;
        N.dst = ccdist;
        neighbors.push_back(N);
      }
    }
  }
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " dataSpheres.txt\n";
  }

  readSpheres(argv[1]);

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("Sphere viewer");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Init the visualizer
  center_x = center_y = center_z = 0.0;  // where we look at
  eye_x = eye_y = 0.0;
  eye_z = 1.0;  // from where we look
  up_x = up_z = 0.0;
  up_y = 1.0;  // direction (normalized)

  mouse_mode = NOTHING;
  view_angle = 45.0;
  znear = 0.01;
  zfar = 10.0;

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_COLOR_MATERIAL);

  // Create light components
  GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
  GLfloat diffuseLight[] = {0.8f, 0.8f, 0.8, 1.0f};
  GLfloat specularLight[] = {0.5f, 0.5f, 0.5f, 1.0f};
  GLfloat positionLight0[] = {1000000.0f, 1000000.0f, 1000000.0f, 1.0f};
  GLfloat positionLight1[] = {-1000000.0f, -1000000.0f, -1000000.0f, 1.0f};

  // Assign created components to GL_LIGHT0
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, positionLight0);

  // Assign created components to GL_LIGHT1
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT1, GL_POSITION, positionLight1);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_BLEND);
  // glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_DEPTH_TEST);
  // glDepthFunc(GL_LEQUAL);

  // ==== Enter GLUT event processing cycle
  adjust_clipping_plans();
  fit_view();
  glutMainLoop();
  return 0;
}

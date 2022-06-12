#ifndef SEE_HPP_E29BD15E
#define SEE_HPP_E29BD15E

//#ifdef __APPLE__
//#include <GLUT/glut.h>
//#else
//#include <GL/glut.h>
//#endif

#include <GL/freeglut.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

int main_window;

int width = 800;
int height = 800;
float wh_ratio = (float)width / (float)height;

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int mouse_start[2];
float view_angle;
float znear;
float zfar;
GLfloat max_length;

double eye_x, eye_y, eye_z;
double center_x, center_y, center_z;
double up_x, up_y, up_z;

int show_background = 0;

struct Sphere {
  double x, y, z;
  double r;
  Sphere() : x(0.0), y(0.0), z(0.0), r(0.0) {}
  Sphere(double x_, double y_, double z_,double r_)
      : x(x_), y(y_), z(z_), r(r_) {}
};

std::vector<Sphere> spheres;

// Box
double xmin, xmax;
double ymin, ymax;
double zmin, zmax;

// Data for drawing spheres
#define X .525731112119133606
#define Z .850650808352039932
static GLfloat vdata[12][3] = {{-X, 0.0, Z}, {X, 0.0, Z},   {-X, 0.0, -Z}, {X, 0.0, -Z}, {0.0, Z, X},  {0.0, Z, -X},
                               {0.0, -Z, X}, {0.0, -Z, -X}, {Z, X, 0.0},   {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
#undef X
#undef Y

static GLuint tindices[20][3] = {{0, 4, 1}, {0, 9, 4},  {9, 5, 4},  {4, 5, 8},  {4, 8, 1},  {8, 10, 1}, {8, 3, 10},
                                 {5, 3, 8}, {5, 2, 3},  {2, 7, 3},  {7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6},
                                 {0, 1, 6}, {6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5},  {7, 2, 11}};

// Drawing functions
void drawsphere(int ndiv, float radius);
void drawtri(GLfloat* a, GLfloat* b, GLfloat* c, int div, float r);
void normalize(GLfloat* a);
void drawParticles();
void drawBoundingBox();

void clear_background();

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);

void rotatePoint(double p_x, double p_y, double p_z, double center_x, double center_y, double center_z, double axis_x,
                 double axis_y, double axis_z, double& ret_x, double& ret_y, double& ret_z, double theta);
void adjust_clipping_plans();
void fit_view();
void bounding_box();
bool fileExists(const char* fileName);
void readSpheres(int num);

#endif /* end of include guard: SEE_HPP_E29BD15E */

# Welcome to the GitHub Pages of `toofus`

`toofus` means `TO`ols `OF`ten `US`ued. It is in fact a collection of C++ header files that I personally develop and use for my own numerical calculation tools. 

A brief description of the different headers and the basic principles of their use is given below.

# Algebra

## `vec2.hpp` and `vec3.hpp`

They are easy-to-use classes for 2D and 3D vectors. Here is example of usage:

```c++
#include <iostream>
#include "vec2.hpp"

int main() {
  vec2r a(1.0, 2.0); // vec2r means vec2<double>
  vec2r b(3.0, 4.0);
  std::cout << "a + b = " << a + b << '\n';
  return 0;
}
```

For all other operations, see directly in the files that are more or less self-documented.

```c++
template <typename T> class vec2 {
public:
  T x, y;

  vec2();
  vec2(T X, T Y);
  vec2(const vec2 &v);
  vec2 &operator=(const vec2 &V) ;
  static vec2 unit_x();
  static vec2 unit_y();
  static vec2 one();
  static vec2 zero();
  void reset();
  void set(T X, T Y);
  void set(T val);
  bool isnull(const T tol = 1e-20) const;
  T *c_vec();
  T &operator[](int i);
  const T &operator[](int i) const;
  const T n() const;
  const T t() const;
  vec2 &operator+=(const vec2 &a);
  vec2 &operator-=(const vec2 &a);
  vec2 &operator*=(T k);
  vec2 &operator/=(T k);
  friend vec2 operator+(const vec2 &a, const vec2 &b);
  friend vec2 operator-(const vec2 &a, const vec2 &b);
  friend vec2 operator-(const vec2 &a);
  friend vec2 operator*(const vec2 &a, T k);
  friend vec2 operator*(T k, const vec2 &a);
  friend vec2 operator/(const vec2 &a, T k);
  friend T operator*(const vec2 &a, const vec2 &b);
  friend vec2<T> component_product(const vec2<T> &a, const vec2<T> &b);
  friend vec2<T> component_min(const vec2<T> &a, const vec2<T> &b);
  friend vec2<T> component_max(const vec2<T> &a, const vec2<T> &b);
  friend vec2<T> component_abs(const vec2<T> &a);
  friend T cross(const vec2<T> &a, const vec2<T> &b);
  friend vec2<T> lerp(double t, const vec2<T> &a, const vec2<T> &b);
  friend T norm2(const vec2 &a);
  friend T norm(const vec2 &a);
  T length() const;
  T normalize();
  vec2 normalized() const;
  friend T determinant(const vec2<T> a, const vec2<T> b);
  bool operator==(const vec2<T> &other) const;
  bool operator!=(const vec2<T> &other) const;
  friend std::ostream &operator<<(std::ostream &pStr, const vec2 &pV);
  friend std::istream &operator>>(std::istream &pStr, vec2 &pV);
};

typedef vec2<double> vec2r;
typedef vec2<int> vec2i;
typedef vec2<unsigned int> vec2ui;
typedef vec2<bool> vec2b;

namespace std {
template <class T> struct less<vec2<T>> { bool operator()(const vec2<T> &lhs, const vec2<T> &rhs) const; };
}
```

## `mat4.hpp` and `mat9.hpp`

They are templated classes for 2-by-2 and 3-by-3 matrices.

```c++
template <typename T> class mat4 {
public:
  T xx, xy;
  T yx, yy;
  
  mat4();
  mat4(const T XX, const T XY, const T YX, const T YY);
  mat4(const T M[]);
  mat4(const mat4 &M);
  mat4 &operator=(const mat4 &M);
  
  static mat4 unit();
  static mat4 zero();
  static mat4 one();
  
  void reset();
  void reset(const double val);
  void set_diag(const double XX, const double YY);
  T &operator[](int i);
  const T &operator[](int i) const;
  T *c_mtx();
  mat4 transposed();
  void transpose();
  void eigenvalues(double &v1, double &v2, bool &swapped) const;
  void eigen(mat4 &V, mat4 &D);
  int sym_eigen(mat4 &V, mat4 &D) const;
  bool inverse();
  mat4 get_inverse();
  double det() const;
  void svd(mat4 &U, mat4 &S, mat4 &V) const;
  bool square_root(mat4 &SqR) const;
  mat4 &operator+=(const mat4 &a);
  mat4 &operator-=(const mat4 &a);
  mat4 &operator*=(double k);
  mat4 &operator/=(double k);
  bool operator==(const mat4 &other) const;
  bool operator!=(const mat4 &other) const;
  friend mat4 operator+(const mat4 &a, const mat4 &b);
  friend mat4 operator-(const mat4 &a, const mat4 &b);
  friend mat4 operator-(const mat4 &a);
  friend mat4 operator*(const mat4 &a, double k);
  friend mat4 operator*(double k, const mat4 &a);
  friend mat4 operator/(const mat4 &a, double k);
  friend vec2r operator*(const mat4 &a, const vec2r &b);
  friend vec2r operator*(const vec2r &b, const mat4 &a);
  friend mat4 operator*(const mat4 &a, const mat4 &b);
  friend std::ostream &operator<<(std::ostream &pStr, const mat4 &pV);
  friend std::istream &operator>>(std::istream &pStr, mat4 &M);
};

typedef mat4<double> mat4r;
typedef mat4<float> mat4f;
typedef mat4<int> mat4i;
typedef mat4<unsigned int> mat4ui;
typedef mat4<bool> mat4b;
```

## `quat.hpp`

It is used for quaternions with the objective of using them for 3D rotations.

```c++
class quat {
public:
  vec3r v;
  double s;
  
  quat();
  quat(double X, double Y, double Z, double S);
  quat(const vec3r &V, const double S);
  quat(const quat &Q);
  static quat identity();
  quat &operator=(const quat &Q);
  quat &operator+=(const quat &a);
  quat &operator-=(const quat &a);
  quat &operator*=(double k);
  quat &operator/=(double k);
  friend quat operator*(const quat &q1, const quat &q2);
  vec3r operator*(const vec3r &V) const;
  quat dot(const vec3r &omega);
  quat ddot(const vec3r &omega, const vec3r &domega);
  void conjugate();
  quat get_conjugated() const;
  void reset();
  void set_axis_angle(const vec3r &V, double angle);
  void set(double X, double Y, double Z, double S);
  double get_angle() const;
  double get_Pitch() const;
  double get_Yaw() const;
  double get_Roll() const;
  vec3r get_axis() const;
  void set_from_to(const vec3r &V1, const vec3r &V2);
  void TwistSwingDecomp(const vec3r &V1, quat &twist, quat &swing);
  void SwingTwistDecomp(const vec3r &V1, quat &swing, quat &twist);
  double normalize();
  void randomize(bool seedTime = false);
  void get_rot_matrix(double M[]) const;
  void get_rot_matrix(mat9r &M) const;
  int set_rot_matrix(double m[]);
  mat9<double> rotate_diag_mat(const vec3r &u) const;
  vec3r rotate(const vec3r &u) const;
  vec3r unrotate(const vec3r &u) const;
  bool operator==(const quat &other) const;
  bool operator!=(const quat &other) const;
  friend std::ostream &operator<<(std::ostream &pStr, const quat &Q);
  friend std::istream &operator>>(std::istream &pStr, quat &Q);
};
```

# Bounding volumes and Space partition

## `AABB.hpp`

An Axis Aligned Bounding Box.

```c++
class AABB {
public:
  vec3r min, max;

  AABB();
  explicit AABB(const vec3r &v);
  AABB(const vec3r &v1, const vec3r &v2);
  AABB(const AABB &aabb);
  explicit AABB(const std::vector<vec3r> &cloud);
  AABB &operator=(const AABB &aabb);
  double getRadius() const;
  void set_single(const vec3r &v);
  void add(const vec3r &v);
  void enlarge(double more);
  void enlarge(const vec3r &more);
  void enlarge(const AABB &more);
  void translate(const vec3r &v);
  bool intersect(const AABB &a) const;
  bool intersect(const vec3r &a) const;
  bool intersectX(const AABB &a) const;
  bool intersectY(const AABB &a) const;
  bool intersectZ(const AABB &a) const;
};
```

## `OBB.hpp`

An Oriented Bounding Box.

## `linkCells.hpp`

```c++
class AABB_Cell {
public:
  std::vector<size_t> bodies;      // Contained bodies
  std::vector<AABB_Cell *> pcells; // Surroundind cells (+ current cell)
};

class linkCells {
public:
  AABB box;                                               // Overall surrounding box
  vec3r minSizes;                                         // Wanted minimum size (along x, y and z) of the cells
  std::vector<std::vector<std::vector<AABB_Cell>>> cells; // Cells that hold only free bodies
  AABB_Cell oversized_bodies;                             // A particular cell that hold only 'too big' bodies
  vec3ui N;                                               // Number of cells in each direction (x, y, and z)

  linkCells(AABB &Box, vec3r &CellMinSizes);
  ~linkCells();
  void init();
  void clear();
  void add_body(size_t B, vec3r &pos, vec3r diag);
  void add_body(size_t B, vec3r &pos, AABB &aabb);
};
```

## `quadtree.hpp` and `octree.hpp`

2D and 3D version of the strategy of partitionning the space with a binary tree. It can fastly query for points inside a box or a disk/sphere. The space-box can also be periodic.

## `OBBtree.hpp`

Space partition with a binary tree of OBBs

# Display

## `ColorTable.hpp`

```c++
#include "../ColorTable.hpp"

#include <iostream>

int main(int argc, char const *argv[]) {
  ColorTable ct;

  ct.setTableID(MATLAB_JET);
  ct.setMinMax(1.5, 21.0);

  colorRGBA col;
  for (double v = 0.0; v <= 25.0; v += 0.5) {
    ct.getRGB(v, &col);
    std::cout << "value: " << v << " -> (" << col.r << ", " << col.g << ", " << col.b << ")\n";
  }
  return 0;
}
```

## `glutTools.hpp`

# Computational Geometry tools

## `geoTool.hpp`

## `geoPack2D.hpp` and `geoPack3D.hpp`

Geometric packing of disks or spheres by means of the Poisson Disk Sampling technique.

```c++
int main(int argc, char const *argv[]) {
  GeoPack2D GP(0.15, 0.15, 5000, 0, 20, 0, 20, 0.0, 0);
  GP.seedTime();
  
  GP.appendDisk(10, 10, 5);
  GP.appendDisk(4, 4, 2);
  GP.appendDisk(16, 4, 2);
  GP.appendDisk(16, 16, 2);
  GP.appendDisk(4, 16, 2);  
  
  GP.reActivate();
  GP.execPeriodic();
  
  GP.parameters(0.05, 0.15, 5000, 0, 20, 0, 20, 0.0, 0);
  GP.reActivate();
  GP.execPeriodic();
  
  GP.save("disks.txt", 5);
  
  return 0;
}
```

## `convexHull.hpp`

See  [this reference](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain)

The function `std::vector<vec2r> convexHull(std::vector<vec2r> &P)` returns a list of points on the convex hull in counter-clockwise order.
Note that the last point in the returned vector is the same as the first one.

```c++
#include <iostream>
#include "convexHull.hpp"

int main (int argc, char const *argv[]) {
        std::vector<vec2r> points;
        points.push_back(vec2r(0,0));
        points.push_back(vec2r(1,1));
        points.push_back(vec2r(1,0));
        points.push_back(vec2r(-1,0));

        std::vector<vec2r> ch = convexHull(points);
        for (size_t i = 0 ; i < ch.size() ; i++) {
                std::cout << ch[i] << std::endl;
        }
        return 0;
}
```

# Maths

## `Mth.hpp`

## `histo.hpp`

```c++
int main (int argc, char const *argv[]) {
  std::vector<double> v;
  for (size_t i = 0 ; i < 1000 ; i++) v.push_back(i/1000.);
  histo H = histo::pdfMaxPerBin(v, 200);
  for (size_t i = 0 ; i < H.data.size() ; i++) 
    std::cout << H.data[i].X << ' ' << H.data[i].ProbDensity << ' ' << H.data[i].Width << '\n';
  return 0;
}
```

## `anova.hpp`

```c++
#include "anova.hpp"

int main(int argc, char const *argv[]) {
  std::vector<std::vector<double>> samples;

  std::vector<double> group1;
  group1.push_back(-7.75405);
  group1.push_back(-7.70286);
  group1.push_back(-7.68725);
  group1.push_back(-7.61047);
  group1.push_back(-7.60942);
  samples.push_back(group1);

  std::vector<double> group2;
  group2.push_back(-7.35701);
  group2.push_back(-7.29485);
  group2.push_back(-7.28961);
  group2.push_back(-7.26047);
  samples.push_back(group2);

  anovaResult result;
  anova::test(samples, result);
  anova::print(result);

  return 0;
}
```

## `linreg.hpp`

## `KStest.hpp`

## `fastSort3.hpp`

# Input / Output

## `message.hpp`

## `fileTool.hpp`

## `ioDataTable.hpp`

## `nextToken.hpp`

## `kwParser.hpp`

# Computation tools

## `factory.hpp`

## `ExecChrono.hpp`


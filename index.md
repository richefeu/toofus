# Welcome to the GitHub Pages of toofus

`toofus` means `TO`ols `OF`ten `US`ued. It is in fact a collection of C++ header files that I personally develop and use for my own numerical calculation tools. 

A brief description of the different headers and the basic principles of their use is given below.

## Algebra

#### `vec2.hpp` and `vec3.hpp`

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

#### `mat4.hpp` and `mat9.hpp`

They are templated classes for 2-by-2 and 3-by-3 matrices.

#### `quat.hpp`

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

## Bounding volumes and Space partition

#### `AABB.hpp`

An Axis Aligned Bounding Box.

```c++
class AABB {
public:
  vec3r min, max;

  AABB();
  explicit AABB(const vec3r &v);
  AABB(const vec3r &v1, const vec3r &v2);
  AABB(const AABB &aabb);
  explicit AABB(const std::vector<vec3r> &cloud) : min(cloud[0]), max(cloud[0]);
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

#### `OBB.hpp`

An Oriented Bounding Box.

#### `linkCells.hpp`

#### `quadtree.hpp` and `octree.hpp`

2D and 3D version of the strategy of partitionning the space with a binary tree. It can fastly query for points inside a box or a disk/sphere. The space-box can also be periodic.

#### `OBBtree.hpp`

Space partition with a binary tree of OBBs

## Display

#### `ColorTable.hpp`

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

#### `glutTools.hpp`

## Computational Geometry

#### `geoTool.hpp`

#### `geoPack2D.hpp` and `geoPack3D.hpp`

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

#### `convexHull.hpp`

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

## Maths

#### `Mth.hpp`

#### `histo.hpp`

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

#### `anova.hpp`

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

#### `linreg.hpp`

#### `KStest.hpp`

#### `fastSort3.hpp`

## Input / Output

#### `message.hpp`

#### `fileTool.hpp`

#### `ioDataTable.hpp`

#### `nextToken.hpp`

#### `kwParser.hpp`

## xxx

#### `factory.hpp`

#### `ExecChrono.hpp`


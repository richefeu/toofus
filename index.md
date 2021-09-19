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

## Bounding volumes and Space partition

#### `AABB.hpp`

An Axis Aligned Bounding Box.

#### `OBB.hpp`

An Oriented Bounding Box.

#### `linkCells.hpp`

#### `quadtree.hpp` and `octree.hpp`

2D and 3D version of the strategy of partitionning the space with a binary tree. It can fastly query for points inside a box or a disk/sphere. The space-box can also be periodic.

#### `OBBtree.hpp`

Space partition with a binary tree of OBBs

## Display

#### `ColorTable.hpp`

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


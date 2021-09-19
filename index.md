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

## Maths

#### `Mth.hpp`

#### `Mth.hpp`

#### `linreg.hpp`

#### `KStest.hpp`

#### `fastSort3.hpp`

## Input / Output

#### `message.hpp`

#### `fileTool.hpp`

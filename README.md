`toofus` brings together a collection of tools in the form of standalone header files that can be included in C++ sources without further complication. 

They are **TO**ols **OF**ten **US**ed in my various simulation codes.


# Setting up


Installation is simply a matter of copying the directory `toofus` somewhere on your computer.

A simple solution to get `toofus` in your home directory, is to copy-past the following command in a terminal (git needs to be installed)

```
cd ~ && rm -rf toofus && git clone https://github.com/richefeu/toofus.git && cd -
```

If `toofus` is installed in e.g. `/usr/local/include`,
since this path should be known by the environment, a header file can simply be included that way:

```c++
#include "toofus/vec3.hpp"
```

It is also possible to simply add the path of the folder at compilation:

```sh
g++ ... -I /path/to/toofus ...
```
For example: `g++ ... -I ~/toofus ...` if the recommended installation has been chosen.

In this case, the inclusion in a C++ source file can be reduced to:

```c++
#include "vec3.hpp"
```

## cmake
The simplest way to integrate `toofus` by means of a `CMakeLists.txt` file is to provide the path to the header files. For example:

```cmake
add_directories(~/toofus)
```

Another solution is to fetch it in the build folder:

```cmake
FetchContent_Declare (
	toofus
  GIT_REPOSITORY https://github.com/richefeu/toofus.git
  GIT_TAG        main
)

FetchContent_GetProperties(toofus)
if(NOT toofus_POPULATED)
  message(STATUS "Fetching toofus")
  FetchContent_Populate(toofus)
  include_directories(${toofus_SOURCE_DIR})
endif()
```

# Avaiblable tools

Not all tools are documented here. These tools are so simple that the best documentation is the code itself. However, the following list allows you to know about what exists. 

### Some maths

The include file `Mth.hpp` provided some useful, easy to use math functions.

### Manipulating 2D / 3D vectors, 2x2 / 3x3 matrices

This is done with the files `vec2.hpp`, `vec3.hpp`, `mat4.hpp` and `mat9.hpp`. 

### Managing the orientation of rigid bodies

A quaternion type is defined in the file `quat.hpp`. It is specifically intended to 3D rotations.


### Axis Aligned Bounding Box (AABB)

An `AABB` is simply two positions (`vec3r`): `min` and `max`. The methods of interest are `intersect(const AABB &a)` and `intersect(const vec3r &a)`.

### Oriented Bounding Box (OBB)

The main method of interest is `intersect(const OBB &obb)`.  It is based on the algorithm described in the book *Real-Time Collision Detection*  by Christer Ericson.

### Data Structures

`linkCells`, `OBBtree.hpp`, `octree.hpp`

### Linear Regression and cubic splines

They are in the files `linreg.hpp` and `cubicSpline.hpp`.

### Histograms, pdf and cdf

They are in the file `histo.hpp` 

### Linear Elasticity

`Rigidity.hpp`, `Compliance.hpp`

# Some other useful header-only libraries


https://github.com/delfrrr/delaunator-cpp


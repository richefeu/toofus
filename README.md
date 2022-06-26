`toofus` brings together a collection of tools in the form of standalone header files that can be included in C++ sources without further complication. 

These are **TO**ols **OF**ten **US**ed in my various simulation codes.

## Setting up


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
For example: `g++ ... -I ~ ...` if the recommended installation has been chosen.

In this case, the inclusion in a C++ source file can be reduced to:

```c++
#include "vec3.hpp"
```


## useful header-only libraries


https://github.com/delfrrr/delaunator-cpp


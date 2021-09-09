`toofus` brings together a collection of tools in the form of standalone header files that can be included in C++ sources without further complication. 

These are **TO**ols **OF**ten **US**ed in my various simulation codes.

## Setting up


Installation is simply a matter of copying the directory `toofus` somewhere on your computer.

Since the path `/usr/local/include` should be known by the environnement, a header file can be included that way:

```c++
#include "toofus/vec3.hpp"
```

It is also possible to simply add the path of the folder at compilation:

```sh
g++ ... -I /path/to/toofus ...
```

In this case, the inclusion can be reduced to:

```c++
#include "vec3.hpp"
```

The recommended and easiest method of installation is to move to the downloaded `toofus` directory containing the header files and run the script `install.sh`:

```sh
sh install.sh
```
## Welcome to the GitHub Pages of toofus

`toofus` means `TO`ols `OF`ten `US`ued. It is in fact a collection of C++ header files that I personally develop and use for my own numerical calculation tools. 

A brief description of the different headers and the basic principles of their use is given below.

### vec2.hpp and vec3.hpp

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

For all other operations, see directly in the file which is self-documented.

### quat.hpp



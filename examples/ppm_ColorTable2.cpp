#include "../ColorTable.hpp"
#include "../ppm.hpp"

#include <iostream>

int main(int argc, char const *argv[]) {

  ColorTable ct;
  colorRGBA col;
  ppm colGrad(20, 250);

  ct.setMinMax(0.4, 24.6);
  ct.setSize(64);

  ct.rebuild_interp_hsv({0, 32, 63}, {{0, 0, 255, 255}, {0, 255, 0, 255}, {255, 0, 0, 255}});

  for (double v = 0.0; v <= 25.0; v += 0.1) {
    ct.getRGB(v, &col);

    for (size_t i = (size_t)(200 * v); i <= (size_t)(200 * v) + 20; i++) {
      colGrad.r[i] = col.r;
      colGrad.g[i] = col.g;
      colGrad.b[i] = col.b;
    }
  }

  colGrad.write("custom.ppm");

  return 0;
}
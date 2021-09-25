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
  ct.savePpm("custom.ppm");

  return 0;
}
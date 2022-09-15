#include "../ColorTable.hpp"

#include <iostream>

int main(int argc, char const *argv[]) {

  ColorTable ct;

  ct.setTableID(MATLAB_JET);
  ct.setMinMax(1.5, 21.0);

  /*
  colorRGBA col;
  for (double v = 0.0; v <= 25.0; v += 0.5) {
    ct.getRGB(v, &col);
    std::cout << "value: " << v << " -> (" << col.r << ", " << col.g << ", " << col.b << ")\n";
  }
  */
  
  color4f col;
  for (double v = 0.0; v <= 25.0; v += 0.5) {
    ct.getColor4f(v, &col);
    std::cout << "value: " << v << " -> (" << col.r << ", " << col.g << ", " << col.b << ")\n";
  }

  return 0;
}
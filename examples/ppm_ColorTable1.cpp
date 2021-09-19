#include "../ColorTable.hpp"
#include "../ppm.hpp"

#include <cstdio>
#include <iostream>

int main(int argc, char const *argv[]) {

  ColorTable ct;
  ppm colGrad(20, 250);

  for (int tableID = 0; tableID <= 20; tableID++) {

    ct.setTableID(tableID);
    ct.setMinMax(0.4, 24.6);
    ct.setSize(64);
    ct.Rebuild();

    colorRGBA col;
    for (double v = 0.0; v <= 25.0; v += 0.1) {
      ct.getRGB(v, &col);

      for (size_t i = (size_t)(200 * v); i <= (size_t)(200 * v) + 20; i++) {
        colGrad.r[i] = col.r;
        colGrad.g[i] = col.g;
        colGrad.b[i] = col.b;
      }
    }

    char fname[256];
    sprintf(fname, "table%d.ppm", tableID);
    colGrad.write(fname);
  }

  return 0;
}
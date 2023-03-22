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
    
    char fname[256];
    snprintf(fname,256, "table%d.ppm", tableID);
    ct.savePpm(fname);

  }

  return 0;
}
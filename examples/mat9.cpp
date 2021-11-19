#include "../mat9.hpp"

#include <iostream>

int main(int argc, char const *argv[]) {
  
  mat9r M;
  M.at(2,2) = 3.;
  
  std::cout << M << "\n";
  std::cout << M.at(2,2) << "\n";
  
  return 0;
}
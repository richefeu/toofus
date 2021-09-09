// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#ifndef IODATATABLE_HPP
#define IODATATABLE_HPP

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

// The resulting table is DataTable[ncol][nlin]
template <typename T>
void loadDataTable(const char * filename, std::vector<std::vector<T> > & DataTable) {
  
  std::ifstream ifs(filename);
  std::string lineStr;
  int ncol = 0;
  int nlin = 0;
  
  while (std::getline(ifs, lineStr)) {
    if (lineStr[0] == '#') continue;
    
    std::istringstream iss(lineStr);
    if (ncol == 0) {
      std::string str(lineStr);
      std::string delims = " \t,";
      std::size_t current, previous = 0;
      current = str.find_first_of(delims);
      ncol++;
      while (current != std::string::npos) {
        previous = current + 1;
        current = str.find_first_of(delims, previous);
        ncol++;
      }
      DataTable.resize(ncol);
    }
    
    T value;
    for (int c = 0 ; c < ncol ; c++) {
      iss >> value;
      DataTable[c].push_back(value);
    }
    nlin++;  
  }
  
  std::cout << "ncol = " << ncol << std::endl;
  std::cout << "nlin = " << nlin << std::endl;
}

template <typename T>
void saveDataTable(const char * filename, std::vector<std::vector<T> > & DataTable) {
  if (DataTable.empty()) return;
  
  std::ofstream ofs(filename);
  int ncol = DataTable.size();
  int nlin = DataTable[0].size();
  for (int l = 0; l < nlin; l++) {
    for (int c = 0; c < ncol; c++) {
      ofs << DataTable[c][l] << ' ';
    }
    ofs << '\n';
  }
}

/*
int main() {
  std::cout << "Load Data Table\n";
  
  std::vector<std::vector<double> > D;
  
  loadDataTable("data.txt", D);
  for (size_t c = 0 ; c < D.size() ; c++) {
    for (size_t i = 0 ; i < D[c].size() ; i++) std::cout << D[c][i] << std::endl;
    std::cout << std::endl;
  }
  
  
  return 0;
}
*/

#endif /* end of include guard: IODATATABLE_HPP */

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

#ifndef PPM_H
#define PPM_H

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Process a binary PPM file
class ppm {

  /**
   * Resets the PPM object to its default state.
   *
   * This method sets the width, height, and maximum color value of the PPM object
   * to their default values of 0, 0, and 255 respectively.
   */
  void init() {
    width       = 0;
    height      = 0;
    max_col_val = 255;
  }

  // info about the PPM file (height and width)
  unsigned int nr_lines;
  unsigned int nr_columns;

public:
  // arrays for storing the R,G,B values
  std::vector<unsigned char> r;
  std::vector<unsigned char> g;
  std::vector<unsigned char> b;

  unsigned int height;      ///< height of the PPM image
  unsigned int width;       ///< width of the PPM image
  unsigned int max_col_val; ///< maximum color value in the PPM image
  unsigned int size;        ///< total number of elements (pixels)

  ppm() {
    init();
  }

  /// constructs a PPM object and loads data from the specified file
  ppm(const std::string &fname) {
    init();
    read(fname);
  }

  /**
   * Creates an "empty" PPM image with the given width and height. The R, G, B
   * arrays are filled with zeros.
   *
   * @param _width The desired width of the PPM image.
   * @param _height The desired height of the PPM image.
   */
  ppm(const unsigned int _width, const unsigned int _height) {
    init();
    width      = _width;
    height     = _height;
    nr_lines   = height;
    nr_columns = width;
    size       = width * height;

    // fill r, g and b with 0
    r.resize(size);
    g.resize(size);
    b.resize(size);
  }

  // read the PPM image from fname
  void read(const std::string &fname) {
    std::ifstream inp(fname.c_str(), std::ios::in | std::ios::binary);
    if (inp.is_open()) {
      std::string line;
      std::getline(inp, line);
      if (line != "P6") {
        std::cout << "Error. Unrecognized file format." << std::endl;
        return;
      }
      std::getline(inp, line);
      while (line[0] == '#') { std::getline(inp, line); }
      std::stringstream dimensions(line);

      try {
        dimensions >> width;
        dimensions >> height;
        nr_lines   = height;
        nr_columns = width;
      } catch (std::exception &e) {
        std::cout << "Header file format error. " << e.what() << std::endl;
        return;
      }

      std::getline(inp, line);
      std::stringstream max_val(line);
      try {
        max_val >> max_col_val;
      } catch (std::exception &e) {
        std::cout << "Header file format error. " << e.what() << std::endl;
        return;
      }

      size = width * height;

      r.reserve(size);
      g.reserve(size);
      b.reserve(size);

      char aux;
      for (unsigned int i = 0; i < size; ++i) {
        inp.read(&aux, 1);
        r[i] = (unsigned char)aux;
        inp.read(&aux, 1);
        g[i] = (unsigned char)aux;
        inp.read(&aux, 1);
        b[i] = (unsigned char)aux;
      }
    } else {
      std::cout << "Error. Unable to open " << fname << std::endl;
    }
    inp.close();
  }

  /**
   * Writes the PPM image to a file. The file name is given in fname.
   * The image is written in binary mode.
   */
  void write(const std::string &fname) {
    std::ofstream inp(fname.c_str(), std::ios::out | std::ios::binary);
    if (inp.is_open()) {

      inp << "P6\n";
      inp << width;
      inp << " ";
      inp << height << "\n";
      inp << max_col_val << "\n";

      char aux;
      for (unsigned int i = 0; i < size; ++i) {
        aux = (char)r[i];
        inp.write(&aux, 1);
        aux = (char)g[i];
        inp.write(&aux, 1);
        aux = (char)b[i];
        inp.write(&aux, 1);
      }
    } else {
      std::cout << "Error. Unable to open " << fname << std::endl;
    }
    inp.close();
  }
};

#endif
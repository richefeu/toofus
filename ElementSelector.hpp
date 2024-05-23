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

#ifndef ELEMENT_SELECTOR_HPP
#define ELEMENT_SELECTOR_HPP

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

template <class T> class ElementSelector {
public:
  std::string inputOption;
  size_t ID;
  std::vector<size_t> ID_LIST;

  double BOX_X0;
  double BOX_Y0;
  double BOX_X1;
  double BOX_Y1;

  double GRID_X0;
  double GRID_Y0;
  double GRID_LX;
  double GRID_LY;
  double GRID_TOL;

  // TODO: use less variables
  /*
  double X0;
  double Y0;
  double X1;
  double Y1;
  double LX;
  double LY;
  double TOL;
  double P1;
  double P2;
  double P3;
  */

  std::function<void(T *)> actionForNone;
  std::function<void(T *, size_t id)> actionForId;
  std::function<size_t(T *)> ContainerSize;
  std::function<void(T *, size_t, double &, double &)> getXY;

  ElementSelector() {
    inputOption = "NONE";
    actionForNone = [](T *) {};
    actionForId = [](T *, size_t) {};
    ContainerSize = [](T *) -> size_t { return 0; };
    getXY = [](T *, size_t, double &x, double &y) {
      x = 0.0;
      y = 0.0;
    };
  }

  void execute(T *Obj = nullptr) {

    if (inputOption == "NONE") {
      actionForNone(Obj);
    } else if (inputOption == "ALL") {
      for (size_t p = 0; p < ContainerSize(Obj); p++) {
        actionForId(Obj, p);
      }
    } else if (inputOption == "ID") {
      actionForId(Obj, ID);
    } else if (inputOption == "LIST") {
      for (size_t i = 0; i < ID_LIST.size(); i++) {
        actionForId(Obj, ID_LIST[i]);
      }
    } else if (inputOption == "GRID") {

      for (size_t p = 0; p < ContainerSize(Obj); p++) {
        double x, y;
        getXY(Obj, p, x, y);
        double X0 = std::round((x - GRID_X0) / GRID_LX) * GRID_LX - 0.5 * GRID_TOL;
        double Y0 = std::round((y - GRID_Y0) / GRID_LY) * GRID_LY - 0.5 * GRID_TOL;
        double X1 = X0 + GRID_TOL;
        double Y1 = Y0 + GRID_TOL;
        if ((x >= X0 && x <= X1) || (y >= Y0 && y <= Y1)) {
          actionForId(Obj, p);
        }
      }

    } else if (inputOption == "BOX") {

      for (size_t p = 0; p < ContainerSize(Obj); p++) {
        double x, y;
        getXY(Obj, p, x, y);
        if (x >= BOX_X0 && x <= BOX_X1 && y >= BOX_Y0 && y <= BOX_Y1) {
          actionForId(Obj, p);
        }
      }
    }
  }

  friend std::istream &operator>>(std::istream &is, ElementSelector &ES) {
    is >> ES.inputOption;
    if (ES.inputOption == "NONE" || ES.inputOption == "ALL") {
    } else if (ES.inputOption == "ID") {
      is >> ES.ID;
    } else if (ES.inputOption == "LIST") {
      size_t nb = 0;
      is >> nb;
      ES.ID_LIST.clear();
      for (size_t i = 0; i < nb; i++) {
        size_t id;
        is >> id;
        ES.ID_LIST.push_back(id);
      }
    } else if (ES.inputOption == "GRID") {
      is >> ES.GRID_X0 >> ES.GRID_Y0 >> ES.GRID_LX >> ES.GRID_LY >> ES.GRID_TOL;
    } else if (ES.inputOption == "BOX") {
      is >> ES.BOX_X0 >> ES.BOX_Y0 >> ES.BOX_X1 >> ES.BOX_Y1;

    } else {
      ES.inputOption = "NONE";
    }

    return is;
  }
};

#endif /* end of include guard: ELEMENT_SELECTOR_HPP */

#if 0

#include "Mth.hpp"
#include "vec2.hpp"

#include <sstream>

struct Elem {
  vec2r pos;
  Elem(double x, double y) {
    pos.x = x;
    pos.y = y;
  }
};

struct AnyBox {
  std::vector<Elem> element;
  AnyBox() {
    for (size_t i = 0; i < 1000; i++) {
      element.push_back(Elem(Mth::random(0., 10.), Mth::random(0., 10.)));
    }
  }
};

int main(int argc, char const *argv[]) {
  AnyBox box;

  ElementSelector<AnyBox> selector;
  selector.actionForNone = [](AnyBox *B) { std::cout << "NONE " << B->element.size() << '\n'; };
  selector.actionForId = [](AnyBox *B, size_t p) {
    std::cout << B->element[p].pos << '\n';
  };
  selector.ContainerSize = [](AnyBox *B) -> size_t { return B->element.size(); };
  selector.getXY = [](AnyBox *B, size_t p, double &x, double &y) {
    x = B->element[p].pos.x;
    y = B->element[p].pos.y;
  };
	
  std::stringstream ss;
  // ss << "LIST 2    2 25";
  // ss << "ALL";
  // ss << "BOX 5 5 9 9";
	ss << "GRID 0 0.2 3 1.5 0.5";

  ss >> selector;
  selector.execute(&box);

  return 0;
}

#endif

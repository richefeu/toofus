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

#ifndef TEMPO_HPP
#define TEMPO_HPP

#include <functional>
#include <string>
#include <vector>

/**
double V1, V2;
Tempo<double> myTempo;
myTempo.plug(&V1);
myTempo.plug(&V2);
myTempo.set("Ramp", 0.0, 1.0, 7.0, 10.0);
myTempo.update(1.0/3.0); // both V1 and V2 will be set to 8.0
*/
template <class T> struct Tempo {
  std::vector<T *> addresses;
  std::function<T(double t)> update;

  Tempo() : update(nullptr) {}
  /*
  Tempo(const Tempo &cpy) {
    for (size_t i = 0 ; i < cpy.addresses.size() ; i++) addresses.push_back(cpy.addresses[i]);
    update = cpy.update;
  }
  */
  ~Tempo() {}

  /**
   @brief Add a target to send the value to.
   @details Add a pointer to the list of targets to send the value to.
   @param[in] add pointer to a value of type T
   */
  void plug(T *add) { addresses.push_back(add); }

  /**
   @brief Send a value to all targets.
   @details For each target that has been plugged (i.e. added to the list of targets to send the value to) with `plug`,
   set the value of the target to the given parameter.
   @param[in] value the value to send to all targets
   */
  void send(T value) {
    for (size_t i = 0; i < addresses.size(); i++) {
      T *add = addresses[i];
      if (add != nullptr)
        *add = value;
    }
  }

  /**
   * @brief Configures the update function based on the specified command.
   *
   * @details This function sets the update function for the Tempo object. The behavior of the update
   * function is determined by the command string. If the command is "Range", the update function
   * sends `p3` to all targets if the input time `t` is between `p1` and `p2`, otherwise it sends `p4`.
   * If the command is "Ramp", the update function performs a linear interpolation between `p3` and `p4`
   * for times between `p1` and `p2`, sending the interpolated value to all targets. Before `p1`, it
   * sends `p3`, and after `p2`, it sends `p4`.
   *
   * @param command The command string that determines the update behavior ("Range" or "Ramp").
   * @param p1 The starting point of the time interval.
   * @param p2 The ending point of the time interval.
   * @param p3 The value to send at the start of the interval or for "Range" mode.
   * @param p4 The value to send at the end of the interval or for "Range" mode.
   */
  void set(const std::string command, double p1 = 0, double p2 = 0, T p3 = 0, T p4 = 0) {
    if (command == "Range") {
      //     ^
      //  p3 |      +-------+
      //     |      |       |
      //  p4 |------+       +-----
      //     +------------------------> t
      //            p1     p2
      update = [this, p1, p2, p3, p4](double t) -> T {
        if (t >= p1 && t <= p2) {
          send(p3);
          return p3;
        }
        send(p4);
        return p4;
      };
    } else if (command == "Ramp") {
      //     ^
      //  p4 |        +---------------
      //     |       /
      //  p3 |------+
      //     +------------------------> t
      //           p1  p2
      update = [this, p1, p2, p3, p4](double t) -> T {
        if (t < p1) {
          send(p3);
          return p3;
        }
        if (t < p2) {
          // linear interpolation
          double v = p3 + (p4 - p3) * (t - p1) / (p2 - p1);
          send(v);
          return v;
        }
        send(p4);
        return p4;
      };
    }
  }
};

#endif /* end of include guard: TEMPO_HPP */

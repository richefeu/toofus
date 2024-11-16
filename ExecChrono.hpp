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

#ifndef EXECCHRONO_HPP
#define EXECCHRONO_HPP

#include <cmath>
#include <ctime>
#include <iostream>
#include <string>

#include "message.hpp"

class ExecChrono {
private:
  std::string m_label;
  clock_t m_start, m_end;

public:
  /// Constructor with default label.
  ///
  /// \code
  /// ExecChrono ec;
  /// ec.start();
  /// // do something
  /// ec.stop();
  /// \endcode
  ExecChrono() : m_label("Measured time") { start(); }

  /// Constructor with custom label.
  ///
  /// \code
  /// ExecChrono ec("Time spent in calculation");
  /// ec.start();
  /// // do something
  /// ec.stop();
  /// \endcode
  ExecChrono(const char *label) : m_label(label) { start(); }

  /// Starts the execution timer by recording the current clock time.
  void start() { m_start = clock(); }

  /// Stops the timer and prints the time elapsed since
  /// \ref start() was called to the standard output, labeled
  /// with the label given in the constructor, in human-readable
  /// format (e.g. "3.2 seconds", "1.3 minutes", etc.).
  void stop() {
    m_end = clock();
    // Now print it
    double ts = (double)(m_end - m_start) / CLOCKS_PER_SEC;
    std::cout << m_label << ": " << msg::HumanReadableSeconds(ts) << std::endl;
  }
};

#endif /* end of include guard: EXECCHRONO_HPP */

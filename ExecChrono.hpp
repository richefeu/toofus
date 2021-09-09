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
  ExecChrono() : m_label("Measured time") { start(); }
  ExecChrono(const char *label) : m_label(label) { start(); }
  void start() { m_start = clock(); }
  void stop() {
    m_end = clock();
    // Now print it
    double ts = (double)(m_end - m_start) / CLOCKS_PER_SEC;
    std::cout << m_label << ": " << msg::HumanReadableSeconds(ts) << std::endl;
  }
};


#endif /* end of include guard: EXECCHRONO_HPP */

#if 0
#include <iostream>
int main (int argc, char const *argv[])
{
	ExecChrono MM; // the chrono will start automatically (but we can use 'start' if we want)
	// Doing something...
	for (int i = 0 ; i < 150 ; i++) {
		std::cout << ".";
	}
	std::cout << '\n';
	MM.stop();
	
	return 0;
}
#endif

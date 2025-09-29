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

#ifndef MESSAGE_HPP
#define MESSAGE_HPP

#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

// For debugging, you can display the value of any variable in the console,
// and the name of the variable is also displayed. Eg. myVar = value
#define __SHOW(V) std::cerr << #V " = " << V << std::flush << std::endl

// This macro ouputs a visible mark. Used for debugging purposes
#define __YOP std::cerr << ">>>>>> o__o <<<<<< LINE: " << __LINE__ << std::flush << std::endl

#define __shouldNeverHappen                                                                                            \
  std::cerr << msg::normal() << msg::fg_red() << msg::bold() << msg::underlined() << "WARN" << msg::resetUnderlined()  \
            << msg::fg_lightRed() << ": "                                                                              \
            << "SHOULD NEVER HAPPEN! " << msg::normal() << "see line " << __LINE__ << " of the file " << __FILE__      \
            << " (func: " << __func__ << ")" << std::endl

#define __printCell(WIDTH, WHAT) std::left << std::setw(WIDTH) << std::setfill(' ') << WHAT << ' '
#define __printNamedCell(WIDTH1, LABEL, WIDTH2, VALUE)                                                                 \
  std::right << std::setw(WIDTH1) << std::setfill(' ') << LABEL << ' ' << std::left << std::setw(WIDTH2)               \
             << std::setfill(' ') << VALUE << ' '

/*
// EXAMPLE
double val = 1.234567890123456789;
double val2 = 56789;
std::cout << __printCell(10,"#val = ") << __printCell(13,val) << __printCell(10,"#val2 = ") << __printCell(13,val2) <<
std::endl; std::cout << __printCell(10,"#val = ") << __printCell(13,val) << __printCell(10,"#val2 = ") <<
__printCell(13,val2) << std::endl; std::cout << __printNamedCell(10,"Value1:",13,val) <<
__printNamedCell(10,"Value2:",13,val2) << std::endl;
*/

class msg {
public:
#if defined(__WIN32) || defined(__WIN64) || defined(__WIN32__)
  static char pathSeparator() {
    return (char)92;
  } // Backslash: '\'
#else
  static char pathSeparator() {
    return (char)47;
  } // Slash: '/'
#endif

  /// @see http://misc.flogisoft.com/bash/tip_colors_and_formatting
  static std::string normal() {
    return std::string("\033[0m");
  }
  static std::string warn() {
    return std::string("\033[1m\033[31m >\033[4mWARN\033[24m: ");
  }
  static std::string info() {
    return std::string("\033[32m >\033[4mINFO\033[24m: ");
  }

  static std::string bold() {
    return std::string("\033[1m");
  }
  static std::string dim() {
    return std::string("\033[2m");
  }
  static std::string underlined() {
    return std::string("\033[4m");
  }
  static std::string blink() {
    return std::string("\033[5m");
  }
  static std::string reverse() {
    return std::string("\033[7m");
  }
  static std::string hidden() {
    return std::string("\033[8m");
  }

  static std::string resetBold() {
    return std::string("\033[21m");
  }
  static std::string resetDim() {
    return std::string("\033[22m");
  }
  static std::string resetUnderlined() {
    return std::string("\033[24m");
  }
  static std::string resetBlink() {
    return std::string("\033[25m");
  }
  static std::string resetReverse() {
    return std::string("\033[27m");
  }
  static std::string resetHidden() {
    return std::string("\033[28m");
  }

  static std::string fg_default() {
    return std::string("\033[39m");
  }
  static std::string fg_black() {
    return std::string("\033[30m");
  }
  static std::string fg_red() {
    return std::string("\033[31m");
  }
  static std::string fg_green() {
    return std::string("\033[32m");
  }
  static std::string fg_yellow() {
    return std::string("\033[33m");
  }
  static std::string fg_blue() {
    return std::string("\033[34m");
  }
  static std::string fg_magenta() {
    return std::string("\033[35m");
  }
  static std::string fg_cyan() {
    return std::string("\033[36m");
  }
  static std::string fg_lightGray() {
    return std::string("\033[37m");
  }
  static std::string fg_darkGray() {
    return std::string("\033[90m");
  }
  static std::string fg_lightRed() {
    return std::string("\033[91m");
  }
  static std::string fg_lightGreen() {
    return std::string("\033[92m");
  }
  static std::string fg_lightYellow() {
    return std::string("\03393[m");
  }
  static std::string fg_lightBlue() {
    return std::string("\033[94m");
  }
  static std::string fg_lightMagenta() {
    return std::string("\033[95m");
  }
  static std::string fg_lightCyan() {
    return std::string("\033[96m");
  }
  static std::string fg_white() {
    return std::string("\033[97m");
  }

  static std::string bg_default() {
    return std::string("\033[49m");
  }
  static std::string bg_black() {
    return std::string("\033[40m");
  }
  static std::string bg_red() {
    return std::string("\033[41m");
  }
  static std::string bg_green() {
    return std::string("\033[42m");
  }
  static std::string bg_yellow() {
    return std::string("\033[43m");
  }
  static std::string bg_blue() {
    return std::string("\033[44m");
  }
  static std::string bg_magenta() {
    return std::string("\033[45m");
  }
  static std::string bg_cyan() {
    return std::string("\033[46m");
  }
  static std::string bg_lightGray() {
    return std::string("\033[47m");
  }
  static std::string bg_darkGray() {
    return std::string("\033[100m");
  }
  static std::string bg_lightRed() {
    return std::string("\033[101m");
  }
  static std::string bg_lightGreen() {
    return std::string("\033[102m");
  }
  static std::string bg_lightYellow() {
    return std::string("\033103[m");
  }
  static std::string bg_lightBlue() {
    return std::string("\033[104m");
  }
  static std::string bg_lightMagenta() {
    return std::string("\033[105m");
  }
  static std::string bg_lightCyan() {
    return std::string("\033[106m");
  }
  static std::string bg_white() {
    return std::string("\033[107m");
  }

  static std::ostream &info(const std::string &message, std::ostream &os = std::cerr) {
    os << normal() << fg_green() << bold() << underlined() << "INFO" << resetUnderlined() << fg_default() << ": "
       << message << normal() << std::endl;
    return os;
  }

  static std::ostream &alert(const std::string &message, std::ostream &os = std::cerr) {
    os << normal() << fg_red() << bold() << underlined() << "WARN" << resetUnderlined() << fg_default() << ": "
       << message << normal() << std::endl;
    return os;
  }

  static std::ostream &strongAlert(const std::string &message, std::ostream &os = std::cerr) {
    os << normal() << fg_red() << bold() << underlined() << "WARN" << resetUnderlined() << fg_lightRed() << ": "
       << message << normal() << std::endl;
    return os;
  }

  static std::ostream &unknown(const std::string &token, std::ostream &os = std::cerr) {
    os << normal() << fg_red() << bold() << underlined() << "WARN" << resetUnderlined() << fg_default() << ": "
       << bold() << token << normal() << " is unknown!" << std::endl;
    return os;
  }

  static std::ostream &unknown(const std::string &token, const std::string &where, std::ostream &os = std::cerr) {
    os << normal() << fg_red() << bold() << underlined() << "WARN" << resetUnderlined() << fg_default() << ": "
       << bold() << fg_green() << token << fg_default() << " @ " << where << normal() << " is unknown!" << std::endl;
    return os;
  }

  // DEPRECATED
  static void yop() {
    std::cerr << "Yop" << std::endl;
  }

  /// Converts a time in seconds into a human readable string, e.g. 1h 30m 45s.
  static std::string HumanReadableSeconds(double ts) {
    std::stringstream ss;
    double ts_hour = floor(ts / 3600.0);
    double ts_min  = floor((ts - ts_hour * 3600.0) / 60.0);
    double ts_sec  = ts - 3600.0 * ts_hour - 60.0 * ts_min;

    if (ts_hour > 0.0) ss << ts_hour << "h " << ts_min << "m " << floor(ts_sec) << "s";
    else if (ts_min > 0.0) ss << ts_min << "m " << ts_sec << "s";
    else ss << ts_sec << "s";
    return ss.str();
  }

  static std::ostream &bestPrecision(std::ostream &os = std::cerr) {
    os << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    return os;
  }

  static std::ostream &normalPrecision(std::ostream &os = std::cerr) {
    os.unsetf(std::ios_base::floatfield); // because std::defaultfloat seem not yet compatible with all gcc compilers
    os << std::setprecision(6);
    return os;
  }

  // Process has done i out of n rounds,
  // and we want a bar of width w and resolution r.
  // USAGE:
  // size_t progress = 0;
  // for (i = 0; i < n; i++) {
  //   ... DOING THINGS ...
  //   msg::loadbar(i, n);
  // }
  static void loadbar(size_t x, size_t n, size_t w = 50, std::ostream &os = std::cerr) {
    if ((x != n) && (x % (n / 100 + 1) != 0)) return;

    float ratio = (float)x / (float)n;
    size_t c    = static_cast<size_t>(ratio * (float)w);

    os << std::setw(3) << (size_t)(ratio * 100) << "% [";
    for (size_t i = 0; i < c; i++) os << "|";
    for (size_t i = c; i < w; i++) os << " ";
    os << "]\r" << std::flush;
  }
};

#endif /* end of include guard: MESSAGE_HPP */

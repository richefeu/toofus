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

// This is a basic logger that use fmt to log things.
// It is not thread-safe
// You can #define FMT_HEADER_ONLY before including this header
// or alternatively, link with fmtlib

// usage example:
/*
#define FMT_HEADER_ONLY
#include "fmtLogger.hpp"
// g++-13 -std=c++11 -I /usr/local/include test.cpp -o test

int main(int argc, char const *argv[]) {
  Logger::setLevel(LogLevel::trace);

  int a = 43;

  Logger::trace("coucou {}", a);
  Logger::debug("coucou {}", a++);
  Logger::info("coucou {}", a);
  Logger::warn("coucou {} {}", a-3, "bilou");
  Logger::error("coucou {}", a);
  Logger::critical("coucou {}", a);

  return 0;
}
*/

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <fmt/core.h>
#include <iostream>

enum class LogLevel { trace, debug, info, warn, error, critical, off };

class Logger {
private:
  static LogLevel level;
  Logger() = delete;
  Logger(const Logger &) = delete;

public:
  static void setLevel(LogLevel t_level) { level = t_level; }

  template <typename... Args> static void trace(fmt::format_string<Args...> fmt, Args &&...args) {
    if (level > LogLevel::trace)
      return;
    fmt::print("[Trace] ");
    fmt::print(fmt, std::forward<Args>(args)...);
    fmt::print("\n");
  }

  template <typename... Args> static void debug(fmt::format_string<Args...> fmt, Args &&...args) {
    if (level > LogLevel::debug)
      return;
    fmt::print("[\033[96mDebug\033[39m] ");
    fmt::print(fmt, std::forward<Args>(args)...);
    fmt::print("\n");
  }

  template <typename... Args> static void info(fmt::format_string<Args...> fmt, Args &&...args) {
    if (level > LogLevel::info)
      return;
    fmt::print("[\033[92mInfo\033[39m] ");
    fmt::print(fmt, std::forward<Args>(args)...);
    fmt::print("\n");
  }

  template <typename... Args> static void warn(fmt::format_string<Args...> fmt, Args &&...args) {
    if (level > LogLevel::warn)
      return;
    fmt::print("[\033[91mWarn\033[39m] ");
    fmt::print(fmt, std::forward<Args>(args)...);
    fmt::print("\n");
  }

  template <typename... Args> static void error(fmt::format_string<Args...> fmt, Args &&...args) {
    if (level > LogLevel::error)
      return;
    fmt::print("[\033[31mError\033[39m] ");
    fmt::print(fmt, std::forward<Args>(args)...);
    fmt::print("\n");
  }

  template <typename... Args> static void critical(fmt::format_string<Args...> fmt, Args &&...args) {
    if (level > LogLevel::critical)
      return;
    fmt::print("[\033[41mCritical\033[49m] ");
    fmt::print(fmt, std::forward<Args>(args)...);
    fmt::print("\n");
  }
};

LogLevel Logger::level = LogLevel::info;

#endif /* end of include guard: LOGGER_HPP */

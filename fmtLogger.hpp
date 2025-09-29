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

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <iostream>

enum class LogLevel { trace, debug, info, warn, error, critical, off };

class Logger {
private:
  inline static LogLevel level = LogLevel::info;
  inline static bool useColors = true; // Control whether to use colors
  Logger()                     = delete;
  Logger(const Logger &)       = delete;

  template <typename... Args>
  static void log(LogLevel msgLevel, const std::string &levelName, const std::string &colorCode,
                  fmt::format_string<Args...> fmt, Args &&...args) {
    if (level > msgLevel) return;

    if (useColors) {
      fmt::print("{}[{}]{} ", colorCode, levelName, "\033[0m");
    } else {
      fmt::print("[{}] ", levelName);
    }

    fmt::print(fmt, std::forward<Args>(args)...);
    fmt::print("\n");
  }

public:
  static void setLevel(LogLevel t_level) {
    level = t_level;
  }

  static void setUseColors(bool useColorsFlag) {
    useColors = useColorsFlag;
  }

  template <typename... Args> static void trace(fmt::format_string<Args...> fmt, Args &&...args) {
    log(LogLevel::trace, "Trace", "\033[37m", fmt, std::forward<Args>(args)...);
  }

  template <typename... Args> static void debug(fmt::format_string<Args...> fmt, Args &&...args) {
    log(LogLevel::debug, "Debug", "\033[96m", fmt, std::forward<Args>(args)...);
  }

  template <typename... Args> static void info(fmt::format_string<Args...> fmt, Args &&...args) {
    log(LogLevel::info, "Info", "\033[92m", fmt, std::forward<Args>(args)...);
  }

  template <typename... Args> static void warn(fmt::format_string<Args...> fmt, Args &&...args) {
    log(LogLevel::warn, "Warn", "\033[93m", fmt, std::forward<Args>(args)...);
  }

  template <typename... Args> static void error(fmt::format_string<Args...> fmt, Args &&...args) {
    log(LogLevel::error, "Error", "\033[91m", fmt, std::forward<Args>(args)...);
  }

  template <typename... Args> static void critical(fmt::format_string<Args...> fmt, Args &&...args) {
    log(LogLevel::critical, "Critical", "\033[41m", fmt, std::forward<Args>(args)...);
  }
};

#endif /* end of include guard: LOGGER_HPP */

#if 0

int main(int argc, char const *argv[]) {
  Logger::setLevel(LogLevel::trace);

  Logger::setUseColors(true);

  int a = 43;

  Logger::trace("coucou {}", a);
  Logger::debug("coucou {}", a++);
  Logger::info("coucou {}", a);
  Logger::warn("coucou {} {}", a-3, "bilou");
  Logger::error("coucou {}", a);
  Logger::critical("coucou {}", a);

  return 0;
}

#endif

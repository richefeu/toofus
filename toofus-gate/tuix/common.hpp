/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file common.hpp
 * @brief Common utilities and platform-specific code for TUIX
 * @author TUIX Team
 * @version 1.0.0
 */
 
#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <memory>
#include <map>
#include <functional>
#include <regex>
#include <limits>

#ifdef _WIN32
    #include <windows.h>
    #include <io.h>
    #define ISATTY _isatty
    #define FILENO _fileno
#else
    #include <unistd.h>
    #include <sys/ioctl.h>
    #define ISATTY isatty
    #define FILENO fileno
#endif

// Configuration macros
#ifndef TUIX_NO_COLOR
    #define TUIX_COLOR_ENABLED
#endif

#ifndef TUIX_DEFAULT_WIDTH
    #define TUIX_DEFAULT_WIDTH 80
#endif

namespace tuix {

// ============================================================================
// UTILITY FUNCTIONS AND TYPES
// ============================================================================

namespace detail {
    
    // Get terminal width
    inline size_t get_terminal_width() {
        #ifdef _WIN32
            CONSOLE_SCREEN_BUFFER_INFO csbi;
            if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi)) {
                return csbi.srWindow.Right - csbi.srWindow.Left + 1;
            }
        #else
            struct winsize w;
            if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0) {
                return w.ws_col;
            }
        #endif
        return TUIX_DEFAULT_WIDTH;
    }

    // Check if terminal supports colors
    inline bool supports_color() {
        #ifdef TUIX_COLOR_ENABLED
            #ifdef _WIN32
                // Enable ANSI escape sequences on Windows 10+
                HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
                DWORD dwMode = 0;
                if (GetConsoleMode(hOut, &dwMode)) {
                    dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
                    SetConsoleMode(hOut, dwMode);
                }
                return ISATTY(FILENO(stdout));
            #else
                return ISATTY(FILENO(stdout));
            #endif
        #else
            return false;
        #endif
    }

    // String utilities
    inline std::string repeat(const std::string& str, size_t count) {
        std::string result;
        result.reserve(str.length() * count);
        for (size_t i = 0; i < count; ++i) {
            result += str;
        }
        return result;
    }

    inline std::vector<std::string> split_lines(const std::string& text) {
        std::vector<std::string> lines;
        std::istringstream stream(text);
        std::string line;
        while (std::getline(stream, line)) {
            lines.push_back(line);
        }
        return lines;
    }

    inline std::string strip_ansi(const std::string& text) {
        std::regex ansi_regex(R"(\x1b\[[0-9;]*m)");
        return std::regex_replace(text, ansi_regex, "");
    }

    inline size_t display_width(const std::string& text) {
        return strip_ansi(text).length();
    }

    inline std::string truncate(const std::string& text, size_t max_width, const std::string& suffix = "...") {
        if (display_width(text) <= max_width) return text;
        
        std::string clean = strip_ansi(text);
        if (clean.length() <= max_width) return text;
        
        if (max_width <= suffix.length()) return suffix.substr(0, max_width);
        return clean.substr(0, max_width - suffix.length()) + suffix;
    }

    // Initialization class
    class initializer {
    public:
        initializer() {
            // Initialize color support
            #ifdef _WIN32
                // Set console to UTF-8 on Windows
                SetConsoleOutputCP(CP_UTF8);
            #endif
        }
    };
    
    // Ensure initialization happens
    static initializer init_instance;
}

} // namespace tuix
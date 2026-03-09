/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file utils.hpp
 * @brief Utility functions
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"
#include "color.hpp"

namespace tuix {

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * @namespace utils
 * @brief Utility functions for terminal applications
 */
namespace utils {

/**
 * @brief Pause execution until user presses a key
 * @param message Message to display (default: "Press any key to continue...")
 */
inline void pause(const std::string& message = "Press any key to continue...") {
    std::cout << message;
    std::cin.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
}

/**
 * @brief Play a beep sound
 */
inline void beep() {
#ifdef _WIN32
    Beep(800, 200);
#else
    std::cout << '\a' << std::flush;
#endif
}

/**
 * @brief Wrap text to a specified width
 * @param text Text to wrap
 * @param width Maximum line width
 * @return Wrapped text
 */
inline std::string wrap_text(const std::string& text, size_t width) {
    if (width == 0 || text.empty()) {
        return text;
    }
    
    std::vector<std::string> lines = detail::split_lines(text);
    std::ostringstream result;
    
    for (size_t i = 0; i < lines.size(); ++i) {
        std::string line = lines[i];
        size_t line_length = detail::display_width(line);
        
        if (line_length <= width) {
            result << line;
        } else {
            size_t start = 0;
            size_t end = 0;
            size_t current_width = 0;
            
            while (end < line.length()) {
                // Skip ANSI escape sequences
                if (line[end] == '\033') {
                    size_t escape_end = line.find('m', end);
                    if (escape_end != std::string::npos) {
                        end = escape_end + 1;
                        continue;
                    }
                }
                
                // Count character width
                size_t char_width = 1; 
                current_width += char_width;
                
                // Check if we need to wrap
                if (current_width > width) {
                    // Find last space before width limit
                    size_t last_space = line.rfind(' ', end - 1);
                    if (last_space != std::string::npos && last_space > start) {
                        result << line.substr(start, last_space - start) << '\n';
                        start = last_space + 1;
                    } else {
                        // No space found, hard break
                        result << line.substr(start, end - start) << '\n';
                        start = end;
                    }
                    current_width = 0;
                }
                
                ++end;
            }
            
            // Add remaining text
            if (start < line.length()) {
                result << line.substr(start);
            }
        }
        
        
        if (i < lines.size() - 1) {
            result << '\n';
        }
    }
    
    return result.str();
}

/**
 * @brief Print a banner with the application name
 * @param app_name Application name
 * @param version Version string
 * @param os Output stream
 */
inline void print_banner(const std::string& app_name, 
                        const std::string& version = "", 
                        std::ostream& os = std::cout) {
    size_t width = detail::get_terminal_width();
    std::string line(width, '=');
    
    os << '\n' << line << '\n';
    
    std::string title = app_name;
    if (!version.empty()) {
        title += " v" + version;
    }
    
    size_t padding = (width - title.length()) / 2;
    os << std::string(padding, ' ') 
       << color::formatter().fg_color(color::fg::bright_blue).format(title) 
       << '\n';
    
    os << line << "\n\n";
}

} // namespace utils

} // namespace tuix
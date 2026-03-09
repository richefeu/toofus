/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file color.hpp
 * @brief Color formatting and styling functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"

namespace tuix {

// ============================================================================
// COLOR SYSTEM
// ============================================================================

namespace color {
    /**
     * @enum fg
     * @brief Foreground color codes for terminal text
     */
    enum class fg : int {
        black = 30, red = 31, green = 32, yellow = 33,
        blue = 34, magenta = 35, cyan = 36, white = 37,
        bright_black = 90, bright_red = 91, bright_green = 92,
        bright_yellow = 93, bright_blue = 94, bright_magenta = 95,
        bright_cyan = 96, bright_white = 97
    };

    /**
     * @enum bg
     * @brief Background color codes for terminal text
     */
    enum class bg : int {
        black = 40, red = 41, green = 42, yellow = 43,
        blue = 44, magenta = 45, cyan = 46, white = 47,
        bright_black = 100, bright_red = 101, bright_green = 102,
        bright_yellow = 103, bright_blue = 104, bright_magenta = 105,
        bright_cyan = 106, bright_white = 107
    };

    /**
     * @enum style
     * @brief Text style codes for terminal text
     */
    enum class style : int {
        reset = 0, bold = 1, dim = 2, italic = 3, underline = 4, blink = 5, reverse = 7, strikethrough = 9
    };

    /**
     * @class formatter
     * @brief Formats text with colors and styles
     */
    class formatter {
    private:
        std::string codes_;
        static bool color_enabled_;

    public:
        formatter() = default;
        
        /**
         * @brief Set the foreground color
         * @param color The color to set
         * @return Reference to this formatter for chaining
         */
        formatter& fg_color(fg color) {
            if (color_enabled_) codes_ += "\033[" + std::to_string(static_cast<int>(color)) + "m";
            return *this;
        }
        
        /**
         * @brief Set the background color
         * @param color The color to set
         * @return Reference to this formatter for chaining
         */
        formatter& bg_color(bg color) {
            if (color_enabled_) codes_ += "\033[" + std::to_string(static_cast<int>(color)) + "m";
            return *this;
        }
        
        /**
         * @brief Add a text style
         * @param s The style to add
         * @return Reference to this formatter for chaining
         */
        formatter& add_style(style s) {
            if (color_enabled_) codes_ += "\033[" + std::to_string(static_cast<int>(s)) + "m";
            return *this;
        }
        
        /**
         * @brief Format text with the configured colors and styles
         * @param text The text to format
         * @return Formatted text string
         */
        std::string format(const std::string& text) const {
            if (!color_enabled_ || codes_.empty()) return text;
            return codes_ + text + "\033[0m";
        }
        
        /**
         * @brief Enable or disable color output
         * @param enabled Whether colors should be enabled
         */
        static void enable_color(bool enabled = true) {
            color_enabled_ = enabled && detail::supports_color();
        }
        
        /**
         * @brief Check if color output is enabled
         * @return True if colors are enabled
         */
        static bool is_enabled() { return color_enabled_; }
    };

    // Static member definition
    bool formatter::color_enabled_ = detail::supports_color();

    // Predefined color functions
    
    /**
     * @brief Format text in red
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string red(const std::string& text) { return formatter().fg_color(fg::red).format(text); }
    
    /**
     * @brief Format text in green
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string green(const std::string& text) { return formatter().fg_color(fg::green).format(text); }
    
    /**
     * @brief Format text in yellow
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string yellow(const std::string& text) { return formatter().fg_color(fg::yellow).format(text); }
    
    /**
     * @brief Format text in blue
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string blue(const std::string& text) { return formatter().fg_color(fg::blue).format(text); }
    
    /**
     * @brief Format text in magenta
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string magenta(const std::string& text) { return formatter().fg_color(fg::magenta).format(text); }
    
    /**
     * @brief Format text in cyan
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string cyan(const std::string& text) { return formatter().fg_color(fg::cyan).format(text); }
    
    /**
     * @brief Format text in white
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string white(const std::string& text) { return formatter().fg_color(fg::white).format(text); }
    
    /**
     * @brief Format text as bold
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string bold(const std::string& text) { return formatter().add_style(style::bold).format(text); }
    
    /**
     * @brief Format text as italic
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string italic(const std::string& text) { return formatter().add_style(style::italic).format(text); }
    
    /**
     * @brief Format text as underlined
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string underline(const std::string& text) { return formatter().add_style(style::underline).format(text); }
    
    /**
     * @brief Format text as dim
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string dim(const std::string& text) { return formatter().add_style(style::dim).format(text); }
    
    /**
     * @brief Format text as blinking
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string blink(const std::string& text) { return formatter().add_style(style::blink).format(text); }
    
    /**
     * @brief Format text as reversed
     * @param text The text to format
     * @return Formatted text string
     */
    inline std::string reverse(const std::string& text) { return formatter().add_style(style::reverse).format(text); }
}

} // namespace tuix
/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file alignment.hpp
 * @brief Text alignment utilities
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"

namespace tuix {

// ============================================================================
// ALIGNMENT HELPERS
// ============================================================================

/**
 * @enum align
 * @brief Text alignment options
 */
enum class align { left, center, right };

/**
 * @brief Align text within a specified width
 * @param text The text to align
 * @param width The width to align within
 * @param alignment The alignment type (left, center, right)
 * @param fill The character to use for padding
 * @return Aligned text string
 */
inline std::string align_text(const std::string& text, size_t width, align alignment = align::left, char fill = ' ') {
    size_t text_width = detail::display_width(text);
    if (text_width >= width) return text;
    
    size_t padding = width - text_width;
    
    switch (alignment) {
        case align::left:
            return text + std::string(padding, fill);
        case align::right:
            return std::string(padding, fill) + text;
        case align::center: {
            size_t left_pad = padding / 2;
            size_t right_pad = padding - left_pad;
            return std::string(left_pad, fill) + text + std::string(right_pad, fill);
        }
    }
    return text;
}

/**
 * @brief Left-align text within a specified width
 * @param text The text to align
 * @param width The width to align within
 * @param fill The character to use for padding
 * @return Left-aligned text string
 */
inline std::string left(const std::string& text, size_t width, char fill = ' ') {
    return align_text(text, width, align::left, fill);
}

/**
 * @brief Right-align text within a specified width
 * @param text The text to align
 * @param width The width to align within
 * @param fill The character to use for padding
 * @return Right-aligned text string
 */
inline std::string right(const std::string& text, size_t width, char fill = ' ') {
    return align_text(text, width, align::right, fill);
}

/**
 * @brief Center text within a specified width
 * @param text The text to align
 * @param width The width to align within
 * @param fill The character to use for padding
 * @return Center-aligned text string
 */
inline std::string center(const std::string& text, size_t width, char fill = ' ') {
    return align_text(text, width, align::center, fill);
}

// ============================================================================
// TITLES AND HEADINGS
// ============================================================================

/**
 * @brief Print a title with decorative fill characters
 * @param title The title text
 * @param width The width of the title line (0 for terminal width)
 * @param fill The character to use for decoration
 * @param os The output stream
 */
inline void print_title(const std::string& title, size_t width = 0, char fill = '=', std::ostream& os = std::cout) {
    if (width == 0) width = detail::get_terminal_width();
    
    size_t title_len = detail::display_width(title);
    if (title_len + 2 >= width) {
        os << title << '\n';
        return;
    }
    
    size_t padding = (width - title_len - 2) / 2;
    std::string left_fill(padding, fill);
    std::string right_fill(width - title_len - 2 - padding, fill);
    
    os << left_fill << ' ' << title << ' ' << right_fill << '\n';
}

/**
 * @brief Print a heading with a specified level
 * @param heading The heading text
 * @param level The heading level (1-6)
 * @param os The output stream
 */
inline void print_heading(const std::string& heading, int level = 1, std::ostream& os = std::cout) {
    std::string prefix(level, '#');
    os << color::bold(prefix + " " + heading) << '\n';
}

// ============================================================================
// DIVIDERS AND SEPARATORS
// ============================================================================

/**
 * @brief Print a horizontal divider line
 * @param width The width of the divider (0 for terminal width)
 * @param symbol The character to use for the divider
 * @param os The output stream
 */
inline void print_divider(size_t width = 0, char symbol = '-', std::ostream& os = std::cout) {
    if (width == 0) width = detail::get_terminal_width();
    os << std::string(width, symbol) << '\n';
}

/**
 * @brief Print a section divider with a label
 * @param label The label text
 * @param width The width of the divider (0 for terminal width)
 * @param symbol The character to use for the divider
 * @param os The output stream
 */
inline void print_section(const std::string& label, size_t width = 0, char symbol = '-', std::ostream& os = std::cout) {
    if (width == 0) width = detail::get_terminal_width();
    
    size_t label_len = detail::display_width(label);
    if (label_len + 6 >= width) {
        os << symbol << symbol << symbol << ' ' << label << ' ' << symbol << symbol << symbol << '\n';
        return;
    }
    
    size_t side_len = (width - label_len - 2) / 2;
    std::string sides(side_len, symbol);
    
    os << sides << ' ' << label << ' ';
    os << std::string(width - label_len - 2 - side_len, symbol) << '\n';
}

/**
 * @brief Create a section divider string with a label
 * @param label The label text
 * @param width The width of the divider (0 for terminal width)
 * @param symbol The character to use for the divider
 * @return Section divider string
 */
inline std::string section(const std::string& label, size_t width = 0, char symbol = '-') {
    if (width == 0) width = detail::get_terminal_width();
    
    size_t label_len = detail::display_width(label);
    if (label_len + 6 >= width) {
        return std::string(3, symbol) + " " + label + " " + std::string(3, symbol);
    }
    
    size_t side_len = (width - label_len - 2) / 2;
    std::string sides(side_len, symbol);
    
    return sides + " " + label + " " + std::string(width - label_len - 2 - side_len, symbol);
}

} // namespace tuix
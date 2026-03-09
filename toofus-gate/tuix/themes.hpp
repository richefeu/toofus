/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file themes.hpp
 * @brief Theme management for consistent styling
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"
#include "color.hpp"

namespace tuix {

// ============================================================================
// THEMES
// ============================================================================

/**
 * @namespace themes
 * @brief Theme management for consistent styling
 */
namespace themes {

/**
 * @struct theme
 * @brief Defines colors and styles for UI elements
 */
struct theme {
    // Text colors
    color::fg text_color = color::fg::white;
    color::fg title_color = color::fg::bright_blue;
    color::fg heading_color = color::fg::blue;
    color::fg success_color = color::fg::green;
    color::fg warning_color = color::fg::yellow;
    color::fg error_color = color::fg::red;
    color::fg info_color = color::fg::cyan;
    
    // Box colors
    color::fg box_border_color = color::fg::blue;
    color::fg box_title_color = color::fg::bright_blue;
    
    // Table colors
    color::fg table_header_color = color::fg::bright_white;
    color::fg table_border_color = color::fg::blue;
    
    // Progress bar colors
    color::fg progress_bar_color = color::fg::green;
    color::fg progress_text_color = color::fg::white;
    
    // Styles
    color::style title_style = color::style::bold;
    color::style heading_style = color::style::bold;
    color::style emphasis_style = color::style::italic;
    color::style highlight_style = color::style::underline;
};

/**
 * @brief Get the default theme
 * @return Default theme with standard colors
 */
inline theme get_default_theme() {
    return theme{};
}

/**
 * @brief Get a dark theme
 * @return Theme optimized for dark backgrounds
 */
inline theme get_dark_theme() {
    theme dark;
    dark.text_color = color::fg::bright_white;
    dark.title_color = color::fg::bright_cyan;
    dark.heading_color = color::fg::cyan;
    dark.box_border_color = color::fg::cyan;
    dark.box_title_color = color::fg::bright_cyan;
    dark.table_border_color = color::fg::cyan;
    dark.table_header_color = color::fg::bright_cyan;
    return dark;
}

/**
 * @brief Get a minimal theme with less colors
 * @return Minimal theme with limited color usage
 */
inline theme get_minimal_theme() {
    theme minimal;
    minimal.text_color = color::fg::white;
    minimal.title_color = color::fg::white;
    minimal.heading_color = color::fg::white;
    minimal.box_border_color = color::fg::white;
    minimal.box_title_color = color::fg::white;
    minimal.table_border_color = color::fg::white;
    minimal.table_header_color = color::fg::white;
    return minimal;
}

// Global current theme
static theme current_theme = get_default_theme();

/**
 * @brief Set the current global theme
 * @param new_theme The theme to set
 */
inline void set_theme(const theme& new_theme) {
    current_theme = new_theme;
}

} // namespace themes

} // namespace tuix
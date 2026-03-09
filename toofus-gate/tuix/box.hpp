/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file box.hpp
 * @brief Box drawing and bordered text functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"
#include "alignment.hpp"

namespace tuix {

// ============================================================================
// BOXED TEXT
// ============================================================================

/**
 * @enum box_style
 * @brief Box border styles
 */
enum class box_style { single, double_, rounded, thick };

namespace detail {
    /**
     * @struct box_chars
     * @brief Characters used for drawing box borders
     */
    struct box_chars {
        char top_left, top_right, bottom_left, bottom_right;
        char horizontal, vertical;
        char top_tee, bottom_tee, left_tee, right_tee, cross;
    };

    /**
     * @brief Get the box characters for a specific style
     * @param style The box style
     * @return Box characters structure
     */
    inline box_chars get_box_chars(box_style style) {
        switch (style) {
            case box_style::single:
                return {'+', '+', '+', '+', '-', '|', '+', '+', '+', '+', '+'};
            case box_style::double_:
                return {'#', '#', '#', '#', '=', '#', '#', '#', '#', '#', '#'};
            case box_style::rounded:
                return {'.', '.', '`', '\'', '-', '|', '+', '+', '+', '+', '+'};
            case box_style::thick:
                return {'*', '*', '*', '*', '=', '*', '*', '*', '*', '*', '*'};
            default:
                return {'+', '+', '+', '+', '-', '|', '+', '+', '+', '+', '+'};
        }
    }
}

/**
 * @class box
 * @brief Creates and renders boxes with borders and optional titles
 */
class box {
private:
    std::vector<std::string> lines_;
    std::string title_;
    box_style style_;
    size_t width_;
    size_t padding_;
    align alignment_;

public:
    /**
     * @brief Construct a new box
     * @param content The content to display in the box
     * @param width The width of the box (0 for auto-width)
     */
    box(const std::string& content, size_t width = 0) 
        : style_(box_style::single), width_(width), padding_(1), alignment_(align::left) {
        if (width_ == 0) width_ = detail::get_terminal_width() - 4;
        lines_ = detail::split_lines(content);
    }
    
    /**
     * @brief Set the box border style
     * @param s The style to use
     * @return Reference to this box for chaining
     */
    box& set_style(box_style s) { style_ = s; return *this; }
    
    /**
     * @brief Set the box title
     * @param t The title text
     * @return Reference to this box for chaining
     */
    box& set_title(const std::string& t) { title_ = t; return *this; }
    
    /**
     * @brief Set the box width
     * @param w The width in characters
     * @return Reference to this box for chaining
     */
    box& set_width(size_t w) { width_ = w; return *this; }
    
    /**
     * @brief Set the internal padding
     * @param p The padding size in characters
     * @return Reference to this box for chaining
     */
    box& set_padding(size_t p) { padding_ = p; return *this; }
    
    /**
     * @brief Set the text alignment within the box
     * @param a The alignment type
     * @return Reference to this box for chaining
     */
    box& set_alignment(align a) { alignment_ = a; return *this; }
    
    /**
     * @brief Render the box as a string
     * @return Formatted box string
     */
    std::string render() const {
        auto chars = detail::get_box_chars(style_);
        std::ostringstream result;
        
        size_t content_width = width_ - 2; // Account for left/right borders
        
        // Top border
        result << chars.top_left;
        if (!title_.empty()) {
            size_t title_space = content_width - 2; // Space for title
            if (detail::display_width(title_) <= title_space) {
                size_t title_len = detail::display_width(title_);
                size_t left_border = (title_space - title_len) / 2;
                size_t right_border = title_space - title_len - left_border;
                
                result << std::string(left_border, chars.horizontal);
                result << ' ' << title_ << ' ';
                result << std::string(right_border, chars.horizontal);
            } else {
                result << std::string(content_width, chars.horizontal);
            }
        } else {
            result << std::string(content_width, chars.horizontal);
        }
        result << chars.top_right << '\n';
        
        // Empty line if padding
        if (padding_ > 0) {
            result << chars.vertical << std::string(content_width, ' ') << chars.vertical << '\n';
        }
        
        // Content lines
        for (const auto& line : lines_) {
            result << chars.vertical;
            std::string padded_content = std::string(padding_, ' ');
            
            size_t available_width = content_width - 2 * padding_;
            std::string formatted_line = detail::truncate(line, available_width);
            formatted_line = align_text(formatted_line, available_width, alignment_);
            
            result << padded_content << formatted_line << padded_content;
            result << chars.vertical << '\n';
        }
        
        // Empty line if padding
        if (padding_ > 0) {
            result << chars.vertical << std::string(content_width, ' ') << chars.vertical << '\n';
        }
        
        // Bottom border
        result << chars.bottom_left;
        result << std::string(content_width, chars.horizontal);
        result << chars.bottom_right << '\n';
        
        return result.str();
    }
    
    /**
     * @brief Output operator for box
     * @param os The output stream
     * @param b The box to output
     * @return The output stream
     */
    /**
     * @brief Convert box to string
     * @return Formatted box string
     */
    std::string to_string() const {
        return render();
    }
    
    /**
     * @brief Output operator for box
     * @param os The output stream
     * @param b The box to output
     * @return The output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const box& b) {
        return os << b.render();
    }
};

/**
 * @brief Print text in a box
 * @param content The content to display
 * @param style The box style
 * @param title Optional title for the box
 * @param os The output stream
 */
inline void print_boxed(const std::string& content, box_style style = box_style::single, 
                       const std::string& title = "", std::ostream& os = std::cout) {
    box b(content);
    b.set_style(style).set_title(title);
    os << b;
}

} // namespace tuix
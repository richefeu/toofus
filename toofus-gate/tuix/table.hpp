/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file table.hpp
 * @brief Table rendering functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include <iostream>
#include<algorithm>
#include "common.hpp"
#include "alignment.hpp"
#include "box.hpp"

namespace tuix {

// ============================================================================
// TABLE RENDERING
// ============================================================================

/**
 * @class table
 * @brief Creates and renders tables with borders, headers, and aligned columns
 */
class table {
private:
    std::vector<std::vector<std::string>> data_;
    std::vector<std::string> headers_;
    mutable std::vector<align> column_alignments_;
    mutable std::vector<size_t> column_widths_;
    box_style style_;
    bool has_headers_;
    size_t padding_;

public:
    /**
     * @brief Construct a new table
     */
    table() : style_(box_style::single), has_headers_(false), padding_(1) {}
    
    /**
     * @brief Set the table headers
     * @param headers Vector of header strings
     * @return Reference to this table for chaining
     */
    table& set_headers(const std::vector<std::string>& headers) {
        headers_ = headers;
        has_headers_ = true;
        column_alignments_.resize(headers.size(), align::left);
        return *this;
    }
    
    /**
     * @brief Add a row to the table
     * @param row Vector of cell strings
     * @return Reference to this table for chaining
     */
    table& add_row(const std::vector<std::string>& row) {
        data_.push_back(row);
        return *this;
    }
    
    /**
     * @brief Set the alignment for a specific column
     * @param col The column index (0-based)
     * @param alignment The alignment type
     * @return Reference to this table for chaining
     */
    table& set_column_alignment(size_t col, align alignment) {
        if (col < column_alignments_.size()) {
            column_alignments_[col] = alignment;
        }
        return *this;
    }
    
    /**
     * @brief Set the table border style
     * @param s The style to use
     * @return Reference to this table for chaining
     */
    table& set_style(box_style s) { style_ = s; return *this; }
    
    /**
     * @brief Set the cell padding
     * @param p The padding size in characters
     * @return Reference to this table for chaining
     */
    table& set_padding(size_t p) { padding_ = p; return *this; }
    
    /**
     * @brief Render the table as a string
     * @return Formatted table string
     */
    std::string render() const {
        if (data_.empty() && !has_headers_) return "";
        
        auto chars = detail::get_box_chars(style_);
        calculate_column_widths();
        
        std::ostringstream result;
        
        // Top border
        render_border(result, chars, true, false, false);
        
        // Headers
        if (has_headers_) {
            render_row(result, headers_, chars);
            render_border(result, chars, false, true, false);
        }
        
        // Data rows
        for (size_t i = 0; i < data_.size(); ++i) {
            render_row(result, data_[i], chars);
        }
        
        // Bottom border
        render_border(result, chars, false, false, true);
        
        return result.str();
    }
    
    /**
     * @brief Output operator for table
     * @param os The output stream
     * @param t The table to output
     * @return The output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const table& t) {
        return os << t.render();
    }

private:
    /**
     * @brief Calculate the width of each column
     */
   // void calculate_column_widths() const;
    
    /**
     * @brief Render a horizontal border
     */
    void render_border(std::ostringstream& result, const detail::box_chars& chars, 
                       bool is_top, bool is_middle, bool is_bottom) const {
        char left = is_top ? chars.top_left : 
                   is_bottom ? chars.bottom_left : chars.left_tee;
        char right = is_top ? chars.top_right : 
                    is_bottom ? chars.bottom_right : chars.right_tee;
        char junction = is_top ? chars.top_tee : 
                       is_bottom ? chars.bottom_tee : chars.cross;
        
        result << left;
        for (size_t i = 0; i < column_widths_.size(); ++i) {
            result << std::string(column_widths_[i], chars.horizontal);
            if (i < column_widths_.size() - 1) {
                result << junction;
            }
        }
        result << right << '\n';
    }
    
    /**
     * @brief Render a table row
     */
    void render_row(std::ostringstream& result, const std::vector<std::string>& row, 
                   const detail::box_chars& chars) const {
        result << chars.vertical;
        for (size_t i = 0; i < column_widths_.size(); ++i) {
            std::string cell = i < row.size() ? row[i] : "";
            align cell_align = i < column_alignments_.size() ? column_alignments_[i] : align::left;
            
            std::string padded_content = std::string(padding_, ' ');
            size_t content_width = column_widths_[i] - 2 * padding_;
            std::string formatted_cell = align_text(cell, content_width, cell_align);
            
            result << padded_content << formatted_cell << padded_content;
            result << chars.vertical;
        }
        result << '\n';
    }
    
    void calculate_column_widths() const {
        size_t num_cols = has_headers_ ? headers_.size() : 
                         (data_.empty() ? 0 : data_[0].size());
        
        column_widths_.clear();
        column_widths_.resize(num_cols, 0);
        
        // Ensure we have enough alignment settings
        if (column_alignments_.size() < num_cols) {
            column_alignments_.resize(num_cols, align::left);
        }
        
        // Check headers
        if (has_headers_) {
            for (size_t i = 0; i < headers_.size(); ++i) {
                column_widths_[i] = std::max<size_t>(column_widths_[i], detail::display_width(headers_[i]));
            }
        }
        
        // Check data
        for (const auto& row : data_) {
            for (size_t i = 0; i < std::min<size_t>(row.size(), num_cols); ++i) {
                column_widths_[i] = std::max<size_t>(column_widths_[i], detail::display_width(row[i]));
            }
        }
        
        // Add padding
        for (auto& width : column_widths_) {
            width += 2 * padding_;
        }
    }
};

/**
 * @brief Print a table with optional headers
 * @param data The table data (rows of cells)
 * @param headers Optional column headers
 * @param os The output stream
 */
inline void print_table(const std::vector<std::vector<std::string>>& data, 
                       const std::vector<std::string>& headers = {}, 
                       std::ostream& os = std::cout) {
    tuix::table t;
    if (!headers.empty()) {
        t.set_headers(headers);
    }
    for (const auto& row : data) {
        t.add_row(row);
    }
    os << t;
}

} // namespace tuix

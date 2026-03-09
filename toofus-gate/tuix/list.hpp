/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file list.hpp
 * @brief Bullet and nested list functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"

namespace tuix {

// ============================================================================
// LISTS
// ============================================================================

/**
 * @brief Print a simple bullet list
 * @param items Vector of list items
 * @param bullet Bullet character to use
 * @param os The output stream
 */
inline void print_bullet_list(const std::vector<std::string>& items, 
                             const std::string& bullet = "•", 
                             std::ostream& os = std::cout) {
    for (const auto& item : items) {
        os << bullet << ' ' << item << '\n';
    }
}

/**
 * @class nested_list
 * @brief Creates and displays hierarchical nested lists
 */
class nested_list {
public:
    /**
     * @struct item
     * @brief Represents a single item in a nested list
     */
    struct item {
        std::string text;                 ///< Item text
        std::vector<item> children;       ///< Child items
        
        /**
         * @brief Construct a new list item
         * @param text Item text
         */
        item(const std::string& text) : text(text) {}
        
        /**
         * @brief Add a child item
         * @param child_text Text for the child item
         * @return Reference to the added child
         */
        item& add(const std::string& child_text) {
            children.emplace_back(child_text);
            return children.back();
        }
    };

private:
    std::vector<item> items_;
    std::vector<std::string> bullets_;
    size_t indent_size_;

public:
    /**
     * @brief Construct a new nested list
     * @param bullets Vector of bullet characters for different levels
     * @param indent_size Number of spaces for each indentation level
     */
    nested_list(const std::vector<std::string>& bullets = {"•", "◦", "▪", "▫"}, size_t indent_size = 2)
        : bullets_(bullets), indent_size_(indent_size) {}
    
    /**
     * @brief Add a top-level item
     * @param text Item text
     * @return Reference to the added item
     */
    item& add(const std::string& text) {
        items_.emplace_back(text);
        return items_.back();
    }
    
    /**
     * @brief Print the list to an output stream
     * @param os The output stream
     */
    void print(std::ostream& os = std::cout) const {
        for (const auto& item : items_) {
            print_item(item, 0, os);
        }
    }

private:
    /**
     * @brief Print a single item and its children
     * @param item The item to print
     * @param level The nesting level
     * @param os The output stream
     */
    void print_item(const item& item, size_t level, std::ostream& os) const {
        std::string indent(level * indent_size_, ' ');
        std::string bullet = bullets_[level % bullets_.size()];
        
        os << indent << bullet << ' ' << item.text << '\n';
        
        for (const auto& child : item.children) {
            print_item(child, level + 1, os);
        }
    }
};

} // namespace tuix
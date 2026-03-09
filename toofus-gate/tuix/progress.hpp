/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file progress.hpp
 * @brief Progress bar functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"

namespace tuix {

// ============================================================================
// PROGRESS BARS
// ============================================================================

/**
 * @class progress_bar
 * @brief Displays and updates a progress bar in the terminal
 */
class progress_bar {
private:
    size_t width_;
    char fill_char_;
    char empty_char_;
    bool show_percentage_;
    std::string prefix_;
    std::string suffix_;

public:
    /**
     * @brief Construct a new progress bar
     * @param width Width of the progress bar in characters
     * @param fill Character to use for filled portion
     * @param empty Character to use for empty portion
     */
    progress_bar(size_t width = 40, char fill = '#', char empty = ' ') 
        : width_(width), fill_char_(fill), empty_char_(empty), show_percentage_(true) {}
    
    /**
     * @brief Set the characters used for the progress bar
     * @param fill Character for filled portion
     * @param empty Character for empty portion
     * @return Reference to this progress bar for chaining
     */
    progress_bar& set_chars(char fill, char empty) {
        fill_char_ = fill; empty_char_ = empty; return *this;
    }
    
    /**
     * @brief Set the prefix text
     * @param prefix Text to display before the progress bar
     * @return Reference to this progress bar for chaining
     */
    progress_bar& set_prefix(const std::string& prefix) { prefix_ = prefix; return *this; }
    
    /**
     * @brief Set the suffix text
     * @param suffix Text to display after the progress bar
     * @return Reference to this progress bar for chaining
     */
    progress_bar& set_suffix(const std::string& suffix) { suffix_ = suffix; return *this; }
    
    /**
     * @brief Enable or disable percentage display
     * @param show Whether to show percentage
     * @return Reference to this progress bar for chaining
     */
    progress_bar& show_percent(bool show) { show_percentage_ = show; return *this; }
    
    /**
     * @brief Render the progress bar as a string
     * @param progress Progress value between 0.0 and 1.0
     * @return Formatted progress bar string
     */
    std::string render(double progress) const {
        progress = (std::max)(0.0, (std::min)(1.0, progress));
        
        size_t filled = static_cast<size_t>(progress * width_);
        size_t empty = width_ - filled;
        
        std::ostringstream result;
        result << prefix_;
        result << '[' << std::string(filled, fill_char_) << std::string(empty, empty_char_) << ']';
        
        if (show_percentage_) {
            result << ' ' << std::setw(3) << static_cast<int>(progress * 100) << '%';
        }
        
        result << suffix_;
        return result.str();
    }
    
    /**
     * @brief Print the progress bar to an output stream
     * @param progress Progress value between 0.0 and 1.0
     * @param os The output stream
     */
    void print(double progress, std::ostream& os = std::cout) const {
        os << '\r' << render(progress) << std::flush;
    }
    
    /**
     * @brief Print a completed progress bar and add a newline
     * @param os The output stream
     */
    void finish(std::ostream& os = std::cout) const {
        print(1.0, os);
        os << '\n';
    }
};

/**
 * @brief Print a simple progress bar
 * @param progress Progress value between 0.0 and 1.0
 * @param width Width of the progress bar in characters
 * @param os The output stream
 */
inline void print_progress(double progress, size_t width = 40, std::ostream& os = std::cout) {
    progress_bar pb(width);
    pb.print(progress, os);
}

} // namespace tuix
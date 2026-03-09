/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file spinner.hpp
 * @brief Spinner animation functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"

namespace tuix {

// ============================================================================
// SPINNERS
// ============================================================================

/**
 * @class spinner
 * @brief Displays an animated spinner in the terminal
 */
class spinner {
private:
    std::vector<std::string> frames_;
    size_t current_frame_;
    std::string message_;
    std::chrono::milliseconds delay_;
    bool running_;
    std::unique_ptr<std::thread> thread_;

public:
    /**
     * @brief Construct a new spinner
     * @param message Text to display next to the spinner
     */
    spinner(const std::string& message = "Loading...") 
        : frames_{"|" , "/", "-", "\\"}, current_frame_(0), message_(message), 
          delay_(100), running_(false) {}
    
    /**
     * @brief Set the animation frames
     * @param frames Vector of frame strings
     * @return Reference to this spinner for chaining
     */
    spinner& set_frames(const std::vector<std::string>& frames) {
        frames_ = frames; return *this;
    }
    
    /**
     * @brief Set the message text
     * @param message Text to display next to the spinner
     * @return Reference to this spinner for chaining
     */
    spinner& set_message(const std::string& message) { message_ = message; return *this; }
    
    /**
     * @brief Set the animation delay
     * @param delay Delay between frames in milliseconds
     * @return Reference to this spinner for chaining
     */
    spinner& set_delay(std::chrono::milliseconds delay) { delay_ = delay; return *this; }
    
    /**
     * @brief Start the spinner animation
     * @param os The output stream
     */
    void start(std::ostream& os = std::cout) {
        if (running_) return;
        running_ = true;
        thread_ = std::make_unique<std::thread>([this, &os]() {
            while (running_) {
                os << '\r' << frames_[current_frame_] << ' ' << message_ << std::flush;
                current_frame_ = (current_frame_ + 1) % frames_.size();
                std::this_thread::sleep_for(delay_);
            }
        });
    }
    
    /**
     * @brief Stop the spinner animation
     * @param os The output stream
     */
    void stop(std::ostream& os = std::cout) {
        if (!running_) return;
        running_ = false;
        if (thread_ && thread_->joinable()) {
            thread_->join();
        }
        os << '\r' << std::string(message_.length() + 10, ' ') << '\r';
        thread_.reset();
    }
    
    /**
     * @brief Destructor ensures the spinner is stopped
     */
    ~spinner() {
        stop();
    }
};

} // namespace tuix
/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file logging.hpp
 * @brief Logging system with colored output and different log levels
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"
#include "color.hpp"

namespace tuix {

// ============================================================================
// LOGGING SYSTEM
// ============================================================================

/**
 * @enum log_level
 * @brief Log severity levels
 */
enum class log_level { debug = 0, info = 1, warning = 2, error = 3, success = 4 };

/**
 * @class logger
 * @brief Logging system with colored output and different log levels
 */
class logger {
private:
    log_level min_level_;
    bool show_timestamp_;
    std::ostream* output_;

public:
    /**
     * @brief Construct a new logger
     * @param min_level Minimum log level to display
     * @param show_timestamp Whether to show timestamps
     * @param output Output stream to write logs to
     */
    logger(log_level min_level = log_level::info, bool show_timestamp = false, std::ostream& output = std::cout)
        : min_level_(min_level), show_timestamp_(show_timestamp), output_(&output) {}
    
    /**
     * @brief Set the minimum log level
     * @param level The minimum level to display
     */
    void set_level(log_level level) { min_level_ = level; }
    
    /**
     * @brief Enable or disable timestamps
     * @param enabled Whether timestamps should be shown
     */
    void set_timestamp(bool enabled) { show_timestamp_ = enabled; }
    
    /**
     * @brief Log a message with a specific level
     * @param level The log level
     * @param message The message to log
     */
    void log(log_level level, const std::string& message) {
        if (level < min_level_) return;
        
        std::string prefix = get_level_prefix(level);
        std::string timestamp = show_timestamp_ ? get_timestamp() : "";
        
        *output_ << timestamp << prefix << message << '\n';
    }
    
    /**
     * @brief Log a debug message
     * @param message The message to log
     */
    void debug(const std::string& message) { log(log_level::debug, message); }
    
    /**
     * @brief Log an info message
     * @param message The message to log
     */
    void info(const std::string& message) { log(log_level::info, message); }
    
    /**
     * @brief Log a warning message
     * @param message The message to log
     */
    void warning(const std::string& message) { log(log_level::warning, message); }
    
    /**
     * @brief Log an error message
     * @param message The message to log
     */
    void error(const std::string& message) { log(log_level::error, message); }
    
    /**
     * @brief Log a success message
     * @param message The message to log
     */
    void success(const std::string& message) { log(log_level::success, message); }

private:
    /**
     * @brief Get the prefix for a log level
     * @param level The log level
     * @return Formatted prefix string
     */
    std::string get_level_prefix(log_level level) const {
        switch (level) {
            case log_level::debug:
                return color::formatter().fg_color(color::fg::cyan).format("[DEBUG] ");
            case log_level::info:
                return color::formatter().fg_color(color::fg::blue).format("[INFO]  ");
            case log_level::warning:
                return color::formatter().fg_color(color::fg::yellow).format("[WARN]  ");
            case log_level::error:
                return color::formatter().fg_color(color::fg::red).format("[ERROR] ");
            case log_level::success:
                return color::formatter().fg_color(color::fg::green).format("[OK]    ");
        }
        return "[LOG]   ";
    }
    
    /**
     * @brief Get a formatted timestamp
     * @return Timestamp string
     */
    std::string get_timestamp() const {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;
        
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&time_t), "%H:%M:%S");
        oss << '.' << std::setfill('0') << std::setw(3) << ms.count() << " ";
        return oss.str();
    }
};

/**
 * @brief Get the global logger instance
 * @return Reference to the global logger
 */
inline logger& get_logger() {
    static logger instance;
    return instance;
}

// Convenience functions

/**
 * @brief Log a debug message
 * @param message The message to log
 */
inline void log_debug(const std::string& message) { get_logger().debug(message); }

/**
 * @brief Log an info message
 * @param message The message to log
 */
inline void log_info(const std::string& message) { get_logger().info(message); }

/**
 * @brief Log a warning message
 * @param message The message to log
 */
inline void log_warning(const std::string& message) { get_logger().warning(message); }

/**
 * @brief Log an error message
 * @param message The message to log
 */
inline void log_error(const std::string& message) { get_logger().error(message); }

/**
 * @brief Log a success message
 * @param message The message to log
 */
inline void log_success(const std::string& message) { get_logger().success(message); }

} // namespace tuix
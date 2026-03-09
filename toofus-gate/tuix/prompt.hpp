/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file prompt.hpp
 * @brief User input prompts and selection functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"
#include "color.hpp"

namespace tuix {

// ============================================================================
// PROMPTS
// ============================================================================

/**
 * @class prompt
 * @brief Handles user input prompts and selections
 */
class prompt {
public:
    /**
     * @brief Get text input from the user
     * @param message The prompt message
     * @param default_value Default value if user enters nothing
     * @return User input string
     */
    static std::string input(const std::string& message, const std::string& default_value = "") {
        std::string input_prompt = message;
        if (!default_value.empty()) {
            input_prompt += " [" + default_value + "]";
        }
        input_prompt += ": ";
        
        std::cout << input_prompt;
        std::string result;
        std::getline(std::cin, result);
        
        if (result.empty() && !default_value.empty()) {
            return default_value;
        }
        
        return result;
    }
    
    /**
     * @brief Get a yes/no confirmation from the user
     * @param message The prompt message
     * @param default_value Default value if user enters nothing (true=yes, false=no)
     * @return User's choice (true=yes, false=no)
     */
    static bool confirm(const std::string& message, bool default_value = true) {
        std::string yes_no = default_value ? "Y/n" : "y/N";
        std::string prompt_text = message + " [" + yes_no + "]: ";
        
        std::cout << prompt_text;
        std::string result;
        std::getline(std::cin, result);
        
        if (result.empty()) {
            return default_value;
        }
        
        // Convert to lowercase for case-insensitive comparison
        std::transform(result.begin(), result.end(), result.begin(), 
                      [](unsigned char c) { return std::tolower(c); });
        
        return (result == "y" || result == "yes" || result == "true");
    }
    
    /**
     * @brief Let the user select from a list of options
     * @param message The prompt message
     * @param options List of options to choose from
     * @param default_index Default selection index (0-based)
     * @return Selected option index (0-based)
     */
    static size_t select(const std::string& message, 
                        const std::vector<std::string>& options,
                        size_t default_index = 0) {
        if (options.empty()) {
            throw std::invalid_argument("Options list cannot be empty");
        }
        
        if (default_index >= options.size()) {
            default_index = 0;
        }
        
        std::cout << message << "\n";
        
        for (size_t i = 0; i < options.size(); ++i) {
            std::cout << "  " << (i + 1) << ". " << options[i] << "\n";
        }
        
        std::string prompt_text = "Enter selection [" + std::to_string(default_index + 1) + "]: ";
        std::cout << prompt_text;
        
        std::string result;
        std::getline(std::cin, result);
        
        if (result.empty()) {
            return default_index;
        }
        
        try {
            size_t selected = std::stoul(result) - 1;
            if (selected < options.size()) {
                return selected;
            }
        } catch (const std::exception&) {
            // Invalid input, return default
        }
        
        return default_index;
    }
};

/**
 * @brief Convenience function to get text input
 * @param message The prompt message
 * @param default_value Default value if user enters nothing
 * @return User input string
 */
inline std::string ask_input(const std::string& message, const std::string& default_value = "") {
    return prompt::input(message, default_value);
}

/**
 * @brief Convenience function to get yes/no confirmation
 * @param message The prompt message
 * @param default_value Default value if user enters nothing (true=yes, false=no)
 * @return User's choice (true=yes, false=no)
 */
inline bool ask_yes_no(const std::string& message, bool default_value = true) {
    return prompt::confirm(message, default_value);
}

} // namespace tuix
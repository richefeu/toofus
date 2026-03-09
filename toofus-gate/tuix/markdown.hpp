/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file markdown.hpp
 * @brief Markdown rendering functionality
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"
#include "color.hpp"

namespace tuix {

// ============================================================================
// MARKDOWN RENDERER
// ============================================================================

/**
 * @class markdown_renderer
 * @brief Renders markdown-like text with formatting
 */
class markdown_renderer {
private:
    std::ostream* output_;

public:
    /**
     * @brief Construct a new markdown renderer
     * @param output Output stream to write to
     */
    markdown_renderer(std::ostream& output = std::cout) : output_(&output) {}
    
    /**
     * @brief Render markdown text
     * @param text The markdown text to render
     */
    void render(const std::string& text) {
        std::istringstream stream(text);
        std::string line;
        bool in_code_block = false;
        
        while (std::getline(stream, line)) {
            // Check for code block markers
            if (line.find("```") == 0) {
                in_code_block = !in_code_block;
                continue;
            }
            
            if (in_code_block) {
                // Code block content
                *output_ << color::formatter().fg_color(color::fg::bright_black).format(line) << '\n';
            } else {
                // Process regular markdown
                process_line(line);
            }
        }
    }
    
    /**
     * @brief Process a single line of markdown
     * @param line The line to process
     */
    void process_line(const std::string& line) {
        // Check for headers
        if (line.size() >= 2 && line[0] == '#') {
            size_t level = 0;
            while (level < line.size() && line[level] == '#') {
                level++;
            }
            
            if (level <= 6 && (level == line.size() || line[level] == ' ')) {
                std::string header_text = (level < line.size()) ? line.substr(level + 1) : "";
                render_header(header_text, level);
                return;
            }
        }
        
        // Process inline formatting
        std::string processed = process_inline_formatting(line);
        *output_ << processed << '\n';
    }
    
    /**
     * @brief Render a header with appropriate formatting
     * @param text Header text
     * @param level Header level (1-6)
     */
    void render_header(const std::string& text, size_t level) {
        auto formatter = color::formatter();
        
        switch (level) {
            case 1:
                formatter.fg_color(color::fg::bright_blue);
                break;
            case 2:
                formatter.fg_color(color::fg::blue);
                break;
            case 3:
                formatter.fg_color(color::fg::bright_white);
                break;
            default:
                formatter.fg_color(color::fg::white);
                break;
        }
        
        *output_ << formatter.format(text) << '\n';
    }
    
    /**
     * @brief Process inline formatting (bold, italic, etc.)
     * @param text Text to process
     * @return Formatted text
     */
    std::string process_inline_formatting(const std::string& text) {
        std::string result = text;
        std::string pattern;
        std::string replacement;
        size_t pos = 0;
        
        // Bold (double asterisks or double underscores)
        pattern = "\\*\\*(.*?)\\*\\*";
        replacement = color::formatter().fg_color(color::fg::bright_white).format("$1");
        result = std::regex_replace(result, std::regex(pattern), replacement);
        
        pattern = "__(.*?)__";
        result = std::regex_replace(result, std::regex(pattern), replacement);
        
        // Italic (single asterisk or single underscore)
        pattern = "\\*(.*?)\\*";
        replacement = color::formatter().fg_color(color::fg::bright_cyan).format("$1");
        result = std::regex_replace(result, std::regex(pattern), replacement);
        
        pattern = "_(.*?)_";
        result = std::regex_replace(result, std::regex(pattern), replacement);
        
        // Underline (between ~ characters)
        pattern = "~(.*?)~";
        replacement = color::formatter().fg_color(color::fg::bright_white).format("$1");
        result = std::regex_replace(result, std::regex(pattern), replacement);
        
        // Inline code (between backticks)
        pattern = "`(.*?)`";
        replacement = color::formatter().fg_color(color::fg::bright_black).format("$1");
        result = std::regex_replace(result, std::regex(pattern), replacement);
        
        return result;
    }
};

/**
 * @brief Stream operator for markdown renderer
 * @param renderer The markdown renderer
 * @param text The markdown text to render
 * @return Reference to the markdown renderer
 */
inline markdown_renderer& operator<<(markdown_renderer& renderer, const std::string& text) {
    renderer.render(text);
    return renderer;
}

/**
 * @brief Print markdown-formatted text
 * @param text The markdown text to print
 * @param os The output stream
 */
inline void print_markdown(const std::string& text, std::ostream& os = std::cout) {
    markdown_renderer renderer(os);
    renderer.render(text);
}

} // namespace tuix
/**
 * TUIX - Terminal User Interface eXtensions
 * 
 * @file terminal.hpp
 * @brief Terminal control functions
 * @author TUIX Team
 * @version 1.0.0
 */

#pragma once

#include "common.hpp"

namespace tuix {

// ============================================================================
// TERMINAL CONTROL
// ============================================================================

/**
 * @namespace terminal
 * @brief Functions for controlling the terminal
 */
namespace terminal {

/**
 * @brief Clear the terminal screen
 */
inline void clear_screen() {
#ifdef _WIN32
    // Windows-specific clear screen
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    DWORD count;
    DWORD cellCount;
    COORD homeCoords = {0, 0};

    if (GetConsoleScreenBufferInfo(hConsole, &csbi)) {
        cellCount = csbi.dwSize.X * csbi.dwSize.Y;
        FillConsoleOutputCharacter(hConsole, ' ', cellCount, homeCoords, &count);
        FillConsoleOutputAttribute(hConsole, csbi.wAttributes, cellCount, homeCoords, &count);
        SetConsoleCursorPosition(hConsole, homeCoords);
    }
#else
    // ANSI escape sequence for clear screen
    std::cout << "\033[2J\033[1;1H";
#endif
}

/**
 * @brief Move the cursor to a specific position
 * @param row Row position (0-based)
 * @param col Column position (0-based)
 */
inline void move_cursor(int row, int col) {
    std::cout << "\033[" << (row + 1) << ";" << (col + 1) << "H";
}

/**
 * @brief Hide the cursor
 */
inline void hide_cursor() {
    std::cout << "\033[?25l";
}

/**
 * @brief Show the cursor
 */
inline void show_cursor() {
    std::cout << "\033[?25h";
}

/**
 * @brief Save the current cursor position
 */
inline void save_cursor() {
    std::cout << "\033[s";
}

/**
 * @brief Restore the previously saved cursor position
 */
inline void restore_cursor() {
    std::cout << "\033[u";
}

/**
 * @brief Get the terminal size
 * @return Pair of (rows, columns)
 */
inline std::pair<int, int> get_size() {
#ifdef _WIN32
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    return {csbi.srWindow.Bottom - csbi.srWindow.Top + 1, 
            csbi.srWindow.Right - csbi.srWindow.Left + 1};
#else
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return {w.ws_row, w.ws_col};
#endif
}

/**
 * @brief Print text at a specific position
 * @param row Row position (0-based)
 * @param col Column position (0-based)
 * @param text Text to print
 */
inline void print_at(int row, int col, const std::string& text) {
    save_cursor();
    move_cursor(row, col);
    std::cout << text;
    restore_cursor();
}

} // namespace terminal

} // namespace tuix
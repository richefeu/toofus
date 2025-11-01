// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

// !!!! c++11 is required !!!!

#ifndef TERMPROMPT_CPP
#define TERMPROMPT_CPP

#include <iostream>
#include <sstream>
#include <string>
#include <termios.h>
#include <unistd.h>
#include <vector>

template <typename T> class Prompt {
private:
  std::string prompt;
  std::string input;
  size_t cursorPos;
  struct termios oldt;

  void setupTerminal() {
    struct termios newt;
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;
    newt.c_lflag &= ~(ICANON | ECHO); // Disable line buffering and echo
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);
  }

  void restoreTerminal() {
    tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
  }

  void moveCursorLeft() {
    if (cursorPos > 0) {
      cursorPos--;
      std::cout << "\033[1D"; // Move cursor left
    }
  }

  void moveCursorRight() {
    if (cursorPos < input.size()) {
      cursorPos++;
      std::cout << "\033[1C"; // Move cursor right
    }
  }

  void insertChar(char c) {
    input.insert(cursorPos, 1, c);
    cursorPos++;
    std::cout << c;
    // Redraw the rest of the input after the cursor
    for (size_t i = cursorPos; i < input.size(); i++) { std::cout << input[i]; }
    // Move cursor back to the correct position
    for (size_t i = cursorPos; i < input.size(); i++) { std::cout << "\033[1D"; }
  }

  void deleteChar() {
    if (cursorPos < input.size()) {
      input.erase(cursorPos, 1);
      // Redraw the rest of the input after the cursor
      std::cout << "\033[0m \033[0m"; // Clear the character
      for (size_t i = cursorPos; i < input.size(); i++) { std::cout << input[i]; }
      std::cout << " \b"; // Clear the last character
      // Move cursor back to the correct position
      for (size_t i = cursorPos; i < input.size(); i++) { std::cout << "\033[1D"; }
    }
  }

  void handleArrowKeys(char c) {
    if (c == 'A') {        // Up arrow (ignored)
    } else if (c == 'B') { // Down arrow (ignored)
    } else if (c == 'C') { // Right arrow
      moveCursorRight();
    } else if (c == 'D') { // Left arrow
      moveCursorLeft();
    }
  }

public:
  Prompt(const std::string &p) : prompt(p), cursorPos(0) {}

  ~Prompt() {
    restoreTerminal();
  }

  T run() {
    setupTerminal();
    // std::cout << prompt << ": ";
    std::cout << "\033[1m" << prompt << "\033[0m" << ": ";
    char c;
    while ((c = getchar()) != '\n') {
      if (c == 127 || c == 8) { // Backspace or Delete
        if (cursorPos > 0) {
          moveCursorLeft();
          deleteChar();
        }
      } else if (c == 27) { // Escape sequence (arrow keys)
        getchar();          // Skip '['
        c = getchar();
        handleArrowKeys(c);
      } else if (c >= 32 && c <= 126) { // Printable ASCII
        insertChar(c);
      }
    }
    std::cout << std::endl;
    restoreTerminal();

    // Parse the input (for numeric types)
    T value;
    std::istringstream iss(input);
    iss >> value;
    return value;
  }
};

#endif
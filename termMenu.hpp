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

#ifndef TERMMENU_CPP
#define TERMMENU_CPP

#include <iostream>
#include <string>
#include <termios.h>
#include <unistd.h>
#include <vector>

class TermMenu {
private:
  std::string prompt;
  std::vector<std::string> choices;
  int selected;
  int startRow;
  int startCol;

public:
  // Default constructor
  TermMenu() : prompt("Please select an option:"), selected(0), startRow(0), startCol(0) {}

  // Constructor that takes a prompt and an initializer list of choices
  TermMenu(const std::string &prompt, std::initializer_list<std::string> choices)
      : prompt(prompt), choices(choices), selected(0), startRow(0), startCol(0) {}

  /// Set the prompt displayed at the top of the menu
  ///
  /// The prompt is displayed with inverted colors, so it
  /// should be short and not too verbose.
  void setPrompt(const std::string &p) {
    prompt = p;
  }

  /// Add a choice to the menu
  ///
  /// The choice is added to the bottom of the menu, and the
  /// selected item is not changed.
  void addChoice(const std::string &choice) {
    choices.push_back(choice);
  }

  /// Display the menu
  ///
  /// This method displays the menu at the top-left of the screen
  /// with the prompt highlighted in bold, and the selected item
  /// highlighted in bold and underlined.
  void displayMenu() {
    std::cout << "\033[" << startRow << ";" << startCol << "H";
    std::cout << "\033[1m" << prompt << "\033[0m" << std::endl;
    for (size_t i = 0; i < choices.size(); ++i) {
      if (i == selected) {
        std::cout << " > \033[1m\033[4m";
      } else {
        std::cout << "   ";
      }
      std::cout << choices[i] << "\033[0m" << std::endl;
    }
  }

  /// Run the menu
  ///
  /// This method runs the menu until the user presses Enter,
  /// and then returns the selected item as an integer.
  ///
  /// The menu is displayed with the prompt highlighted in bold,
  /// and the selected item is highlighted in bold and underlined.
  ///
  /// The menu is displayed in the middle of the screen, or
  /// as close to the middle as possible.
  ///
  /// The menu can be navigated using the arrow keys.
  int run() {
    char c;
    struct termios oldt, newt;
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;
    newt.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);
    int rows, cols;
    std::string cmd = "stty size";
    FILE *pipe      = popen(cmd.c_str(), "r");
    if (pipe == nullptr) {
      std::cerr << "Error: could not run stty command" << std::endl;
      return -1;
    }
    if (fscanf(pipe, "%d %d", &rows, &cols) != 2) {
      std::cerr << "Error: could not parse stty output" << std::endl;
      pclose(pipe);
      return -1;
    }
    pclose(pipe);
    std::cout << "\033[6n";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '[');
    std::cin >> startRow;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), ';');
    std::cin >> startCol;
    if (startRow > rows - choices.size() - 1) {
      startRow = rows - choices.size() - 1;
      for (size_t l = 0; l < choices.size() + 1; l++) { std::cout << '\n'; }
    }
    displayMenu();
    while ((c = getchar()) != '\n') {
      if (c == '\033') {
        getchar();
        c = getchar();
        if (c == 'A') {
          selected = (selected - 1 + choices.size()) % choices.size();
        } else if (c == 'B') {
          selected = (selected + 1) % choices.size();
        }
      }
      displayMenu();
    }
    tcsetattr(STDIN_FILENO, TCSANOW, &oldt);

    return selected + 1;
  }
};

#endif /* end of include guard: TERMMENU_CPP_44369901 */

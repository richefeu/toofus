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

#ifndef CONSOLE_PROGRESS_BAR
#define CONSOLE_PROGRESS_BAR

#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>

class ConsoleProgressBar {
private:
  size_t nmax;
  size_t width;
  std::string title;
  char progressChar;
  char voidChar;
  char openChar;
  char closeChar;

  /**
   * Validates a character to ensure it's within the printable ASCII range.
   * If the character `c` is a printable character (ASCII 33 to 127 inclusive),
   * it is returned as is. Otherwise, the default character `cdef` is returned.
   *
   * @param c The character to validate.
   * @param cdef The default character to return if `c` is not valid.
   * @return The validated character, either `c` or `cdef`.
   */
  char validChar(char c, char cdef) {
    int ci = (int)c;
    if (ci >= 33 && ci <= 127) {
      return (char)c;
    }
    return (char)cdef;
  }

public:

  /**
   * Constructs a ConsoleProgressBar object with default or specified settings.
   *
   * Initializes the progress bar with a maximum value of `n`, a default width of 50,
   * an empty title, and default characters for progress representation.
   *
   * @param n The maximum value of the progress bar (default is 100).
   *          This determines the completion point of the progress bar.
   */
  ConsoleProgressBar(size_t n = 100) {
    nmax = n;
    width = 50;
    title = "";
    progressChar = '|';
    voidChar = ' ';
    openChar = '[';
    closeChar = ']';
  }

  void setMax(size_t n) { nmax = n; }
  void setTitle(const char *t) { title = t; }
  void setWidth(size_t w) { width = w; }

  void setProgressChar(char c) { progressChar = (char)validChar(c, '|'); }
  void setVoidChar(char c) { voidChar = (char)validChar(c, ' '); }
  void setOpenChar(char c) { openChar = (char)validChar(c, '['); }
  void setCloseChar(char c) { closeChar = (char)validChar(c, ']'); }
  
  /**
   * Updates the progress bar by displaying the current status.
   *
   * @param x The current value of the progress bar (should be between 0 and `nmax`).
   * @param os The output stream to use for displaying the progress bar
   *          (default is `std::cerr`).
   *
   * The `update` method displays the current status of the progress bar.
   * It prints the title, the percentage completed, and the progress bar itself.
   * The progress bar is displayed as a sequence of `progressChar` characters
   * followed by `voidChar` characters.
   *
   * If the progress bar is complete, the `update` method will not print anything.
   * Otherwise, it will print the progress bar to the specified output stream.
   */
  void update(size_t x, std::ostream &os = std::cerr) {
    if ((x != nmax) && (x % (nmax / 100 + 1) != 0))
      return;

    float ratio = (float)x / (float)nmax;
    size_t c = (size_t)std::floor(ratio) * width;

    os << title << std::setw(3) << (size_t)(ratio * 100) << "% " << openChar;
    for (size_t i = 0; i < c; i++) {
      os << progressChar;
    }
    for (size_t i = c; i < width; i++) {
      os << voidChar;
    }
    os << closeChar << '\r' << std::flush;
  }
};

#if 0
int main(int argc, char const *argv[]) {
  ConsoleProgressBar cpb(100);
  cpb.setTitle("Title: ");
  cpb.setWidth(50);
  cpb.setProgressChar(124);
  cpb.setVoidChar('-');
  cpb.update(10, std::cout);
  return 0;
}
#endif

#endif /* end of include guard: CONSOLE_PROGRESS_BAR */

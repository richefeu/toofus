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

#ifndef GLUT_TOOLS_HPP
#define GLUT_TOOLS_HPP

#include "glTools.hpp"


// FIXME: remane glutInputBox
class glInputBox {
public:
  bool inputMode;
  size_t cursorPos;
  int *width;
  int *height;

  std::string inputText;
  std::string promptText;
  std::function<void()> onEnterCallback;

  glInputBox(int *W, int *H) : inputMode(false), cursorPos(0), width(W), height(H) {}

  void start(const std::string &prompt, std::function<void()> callback) {
    inputMode = true;
    promptText = prompt;
    inputText.clear();
    cursorPos = 0;
    onEnterCallback = callback;
  }

  void keyboard(unsigned char key) {
    switch (key) {
    case 13: // Enter key
      inputMode = false;
      onEnterCallback();
      break;
    case 8: // Backspace key
      if (cursorPos > 0) {
        inputText.erase(cursorPos - 1, 1);
        --cursorPos;
      }
      break;
    default: {
      inputText.insert(cursorPos, 1, key);
      ++cursorPos;
    } break;
    }
    glutPostRedisplay(); // Redraw the scene to update the display
  }

  void specialKeyboard(int key) {
    switch (key) {
    case GLUT_KEY_LEFT:
      if (cursorPos > 0) {
        --cursorPos;
      }
      break;
    case GLUT_KEY_RIGHT:
      if (cursorPos < inputText.size()) {
        ++cursorPos;
      }
      break;
    }
    glutPostRedisplay();
  }

  void render() {
    switch2D::go(*width, *height);

    glColor4f(0.91f, 0.52f, 0.6f, 1.0f);
    glBegin(GL_QUADS);
    glVertex2i(0, *height);
    glVertex2i(*width, *height);
    glVertex2i(*width, *height - 22);
    glVertex2i(0, *height - 22);
    glEnd();

    glColor3i(0, 0, 0);
    char txt[256];
    snprintf(txt, 256, "%s: %s", promptText.c_str(), inputText.c_str());
    glText::print(4, *height - 16, txt);
    int shiftPrompt = static_cast<int>(promptText.length()) + 2;
    int shift = 4 + (shiftPrompt + static_cast<int>(cursorPos)) * 10;
    glText::print(shift, *height - 16, "_");
    
    switch2D::back();
  }
};


#endif /* end of include guard: GLUT_TOOLS_HPP */

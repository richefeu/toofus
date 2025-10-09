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

#ifndef FILETOOL_HPP
#define FILETOOL_HPP

#include <string>
#include <cctype>

#if defined(__WIN32) || defined(__WIN64) || defined(__WIN32__)
#include <direct.h> // mkdir
#include <io.h>     // access
#endif

#include <fstream>
#include <iostream>

#include <sys/stat.h> // For mkdir (linux)
#include <unistd.h>   // For access

class fileTool {
public:
  
  static char separator() {
#if defined(__WIN32) || defined(__WIN64) || defined(__WIN32__)
    return 92; // Backslash: '\'
#else
    return 47; // Slash: '/'
#endif
  }

  static bool containsOnlyDigits(const char *str) {
    while (*str != '\0') {       // Iterate through each character in the string
      if (!std::isdigit(*str)) { // Check if the current character is not a digit
        return false;            // If any non-digit character is found, return false
      }
      str++; // Move to the next character
    }
    return true; // Return true if all characters are digits
  }



  /// @brief Robust and portable function to test if a file exists.
  ///
  /// @param[in] fileName name of the file to check.
  ///
  /// @return true if the file exists, false otherwise.
  static bool fileExists(const char *fileName) {
    std::fstream fin;
    fin.open(fileName, std::ios::in);
    if (fin.is_open()) {
      fin.close();
      return true;
    }
    fin.close();
    return false;
  }

  static bool fileExists(const std::string &fileName) {
    return fileExists(fileName.c_str());
  }

  /**
   * @brief Extracts the file extension from a given file name.
   *
   * This function searches for the last occurrence of the dot character ('.')
   * in the provided file name string and returns the substring following it,
   * which represents the file extension. If no dot is found, an empty string
   * is returned.
   *
   * @param FileName The name of the file from which to extract the extension.
   * @return The file extension as a string, or an empty string if no extension is found.
   */
  static std::string GetFileExt(const std::string &FileName) {
    std::string::size_type s = FileName.find_last_of('.');
    if (s != std::string::npos) return FileName.substr(s + 1);
    return std::string("");
  }

  /**
   * @brief Extracts the file path from a given file name.
   *
   * This function searches for the last occurrence of the path separator
   * character (either '/' or '\\') in the provided file name string and
   * returns the substring preceding it, which represents the file path. If
   * no separator is found, the current working directory ("." or ".\")
   * is returned.
   *
   * @param FileName The name of the file from which to extract the path.
   * @return The file path as a string, or the current working directory if no path is found.
   */
  static std::string GetFilePath(const std::string &FileName) {
    std::string::size_type s = FileName.find_last_of(separator());
    if (s != std::string::npos) return FileName.substr(0, s);
    return std::string(".");
  }

  /**
   * @brief Extracts the file name without extension from a given file name.
   *
   * This function searches for the last occurrence of the path separator
   * character (either '/' or '\\') and the dot character ('.') in the
   * provided file name string and returns the substring between them,
   * which represents the file name without extension. If no separator is
   * found, the full file name is returned.
   *
   * @param FileName The name of the file from which to extract the file name.
   * @return The file name without extension as a string, or the full file name if no path is found.
   */
  static std::string GetFileName(const std::string &FileName) {
    std::string::size_type s1 = FileName.find_last_of(separator());
    std::string::size_type s2 = FileName.find_last_of('.');
    if (s1 != std::string::npos && s2 != std::string::npos) return FileName.substr(s1 + 1, s2 - s1 - 1);
    else if (s1 == std::string::npos && s2 != std::string::npos) return FileName.substr(0, s2);
    else if (s1 != std::string::npos && s2 == std::string::npos) return FileName.substr(s1 + 1);
    return FileName;
  }

  /**
   * @brief Creates a directory with the given name.
   *
   * This overloaded version of create_folder() takes a std::string argument
   * and calls the const char * version with the c_str() method.
   *
   * @param folder The name of the directory to be created.
   */
  static void create_folder(std::string &folder) {
    create_folder(folder.c_str());
  }

  /**
   * @brief Creates a directory with the given name if it does not exist.
   *
   * This function checks if a directory with the specified name already exists.
   * If it does not, it attempts to create the directory using platform-specific
   * calls to `mkdir`. On Windows, it uses the default `mkdir` function. On UNIX-like
   * systems, it uses `mkdir` with specific permissions allowing the owner to read,
   * write, and execute, while group and others can read and execute.
   *
   * @param folder The name of the directory to be created as a C-style string.
   */
  static void create_folder(const char *folder) {
    if (access(folder, F_OK)) {
      int stat;
#if defined(__WIN32) || defined(__WIN64) || defined(__WIN32__)
      stat = mkdir(folder);
#else
      stat = mkdir(folder, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif
      if (stat == -1) std::cout << "Cannot create the folder " << folder << std::endl;
    }
  }
};

#endif /* end of include guard: FILETOOL_HPP */

#if 0
#include <iostream>

int main(int argc, char const *argv[]) {
  if (fileTool::fileExists("./examples")) {
    std::cout << "yep\n";
  } else {

    std::cout << "not yep\n";
  }
  return 0;
}
#endif

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

  /// @brief Robust and portable function to test if a file exists
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

  static bool fileExists(const std::string &fileName) { return fileExists(fileName.c_str()); }

  static std::string GetFileExt(const std::string &FileName) {
    std::string::size_type s = FileName.find_last_of('.');
    if (s != std::string::npos)
      return FileName.substr(s + 1);
    return std::string("");
  }

  // Remark: without the last '/' or '\'
  static std::string GetFilePath(const std::string &FileName) {
    std::string::size_type s = FileName.find_last_of(separator());
    if (s != std::string::npos)
      return FileName.substr(0, s);
    return std::string(".");
  }

  static std::string GetFileName(const std::string &FileName) {
    std::string::size_type s1 = FileName.find_last_of(separator());
    std::string::size_type s2 = FileName.find_last_of('.');
    if (s1 != std::string::npos && s2 != std::string::npos)
      return FileName.substr(s1 + 1, s2 - s1 - 1);
    else if (s1 == std::string::npos && s2 != std::string::npos)
      return FileName.substr(0, s2);
    else if (s1 != std::string::npos && s2 == std::string::npos)
      return FileName.substr(s1 + 1);
    return FileName;
  }

  static void create_folder(std::string &folder) {
    create_folder(folder.c_str());
  }
  
  static void create_folder(const char *folder) {
    if (access(folder, F_OK)) {
      int stat;
#if defined(__WIN32) || defined(__WIN64) || defined(__WIN32__)
      stat = mkdir(folder);
#else
      stat = mkdir(folder, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif
      if (stat == -1)
        std::cout << "Cannot create the folder " << folder << std::endl;
    }
  }
  
};

#endif /* end of include guard: FILETOOL_HPP */

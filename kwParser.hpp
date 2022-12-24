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

#ifndef KWPARSER_CPP
#define KWPARSER_CPP

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <sstream>

#define __DO__(IS) [&](std::istream & IS)
#define __GET__(IS, WHAT) [&](std::istream &IS) { IS >> (WHAT); }
class kwParser {
public:
  kwParser() : warn(true) { breakStr = "EOF"; }
  kwParser(bool Warn) : warn(Warn) { breakStr = "EOF"; }

  // This function is actually not very usefull
  // instead do:
  //    myParser.kwMap["key"] = [&](std::istream & is) { is >> values; };
  // or (simpler)
  //    myParser.kwMap["key"] = __DO__ { is >> values; };
  void addKw(std::string kw, std::function<void(std::istream &)> func) { kwMap[kw] = func; }

  void parse(const char *filename) {
    std::ifstream file(filename);
    if (file)
      parse(file);
    else
      std::cerr << "@kwParser::parse, file " << filename << " cannot be opened" << std::endl;
  }
  
  void parseString(const char *chain) {
    std::stringstream ss;
    ss.str(chain);
    if (ss)
      parse(ss);
    else
      std::cerr << "@kwParser::parseString, stream cannot be created" << std::endl;
  }

  void parse(std::istream &is) {
    std::string token;
    is >> token;
    while (is.good()) {
      if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
        getline(is, token); // ignore the rest of the current line
        is >> token;        // next token
        continue;
      } else if (token == breakStr)
        break;
      else {
        auto it = kwMap.find(token);
        if (it != kwMap.end())
          it->second(is);
        else if (warn) {
          std::cerr << "@kwParser::parse, unknown token: " << token << std::endl;
        }
      }

      is >> token;
    }
  }

  std::map<std::string, std::function<void(std::istream &)>> kwMap;
  std::string breakStr;
  bool warn;
};

// EXAMPLE OF USAGE

/**
@file

Text file to be parsed:
@code
# Here are the data to be parsed
myClass.value 123.456
myClass.value2 654.321
myClass.str coucou
# myClass.str coucouComment
@endcode
*/

#if 0
//Example of usage:
struct myClass {
        double value;
        double value2;
        std::string str;
} mc;

int main ()
{
        kwParser parser(false);
        parser.kwMap["myClass.value"] = __DO__(is) { is >> mc.value; };
        parser.kwMap["myClass.value2"] = __GET__(is, mc.value2);
        parser.addKw("myClass.str",[&](std::istream & is) { is >> mc.str; } );

        parser.parse("kwParser.hpp");

        std::cout << "myClass.value = " << mc.value << std::endl;
        std::cout << "myClass.value2 = " << mc.value2 << std::endl;
        std::cout << "myClass.str = " << mc.str << std::endl;
        
        parser.parseString("myClass.str toto");
        std::cout << "myClass.str = " << mc.str << std::endl;
        
        parser.parseString("myClass.value2 6.1 myClass.value 13.6");
        std::cout << "myClass.value = " << mc.value << std::endl;
        std::cout << "myClass.value2 = " << mc.value2 << std::endl;

        return 0;
}
#endif

#endif /* end of include guard: KWPARSER_CPP */

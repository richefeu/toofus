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

/*
 * ConfigJson Tool
 *
 * Overview:
 * The ConfigJson tool is a lightweight and straightforward utility designed for parsing and manipulating JSON data in
 * C++. It provides basic functionality to read, write, and access JSON data, making it simple to integrate into
 * projects that require minimal JSON handling capabilities.
 *
 * Key Features:
 * - Simplicity: The tool is designed to be easy to use, with a simple API for parsing JSON files and accessing data.
 *   It allows users to quickly read JSON data into a structured format and retrieve values using intuitive methods.
 * - Basic JSON Support: It supports fundamental JSON structures, including objects and arrays, and can handle basic
 *   data types such as strings and numbers.
 * - No External Dependencies: The tool is self-contained, relying only on standard C++ libraries, which simplifies
 *   integration and deployment.
 *
 * Limitations:
 * - Limited Data Types: Unlike full-featured JSON parsers, this tool primarily focuses on strings and numbers.
 *   It lacks support for more complex data types such as booleans and null values.
 * - Error Handling: The tool uses basic error handling mechanisms, such as console output for errors, rather than
 *   sophisticated error recovery or detailed error reporting.
 * - Performance: It may not be optimized for performance-critical applications. Parsing large JSON files could be
 *   slower compared to dedicated libraries optimized for speed and efficiency.
 * - No Schema Validation: The tool does not support JSON schema validation, which is essential for applications
 *   requiring strict data structure compliance.
 * - Limited Functionality: Advanced features such as JSON path queries, complex manipulations, or transformations
 *   are not supported.
 * - Basic Parsing: The parsing logic is relatively simplistic and may not handle all edge cases or malformed JSON
 *   data gracefully.
 *
 * Usage:
 * This tool is suitable for simple applications or projects where ease of use and quick setup are prioritized over
 * advanced features and performance. For more complex applications requiring robust JSON handling, a dedicated JSON
 * library would be more appropriate.
 */

/*
 * ConfigJson Tool Usage Documentation
 *
 * Overview:
 * The ConfigJson tool is designed to provide a simple interface for parsing and manipulating JSON data in C++.
 * It is ideal for applications requiring basic JSON handling without the complexity of a full-featured JSON library.
 *
 * Usage Examples:
 *
 * Example 1: Parsing a JSON File
 *
 * ```cpp
 * #include "config_json.hpp"
 * #include <iostream>
 *
 * int main() {
 *     ConfigJson config;
 *     config.parseJson("example.json"); // Parse a JSON file
 *
 *     // Accessing values
 *     std::cout << "Name: " << std::get<std::string>(config["name"]) << std::endl;
 *     std::cout << "Age: " << static_cast<int>(std::get<jsonNumber>(config["age"])) << std::endl;
 *
 *     return 0;
 * }
 * ```
 *
 * Example 2: Modifying and Saving JSON Data
 *
 * ```cpp
 * #include "config_json.hpp"
 * #include <iostream>
 *
 * int main() {
 *     ConfigJson config;
 *     config.parseJson("example.json");
 *
 *     // Modify a value
 *     config["age"] = jsonNumber(30);
 *
 *     // Save the modified configuration to a new file
 *     config.save("output.json");
 *     std::cout << "Configuration saved to output.json" << std::endl;
 *
 *     return 0;
 * }
 * ```
 *
 * Example 3: Accessing Nested JSON Objects
 *
 * ```cpp
 * #include "config_json.hpp"
 * #include <iostream>
 *
 * int main() {
 *     ConfigJson config;
 *     config.parseJson("example.json");
 *
 *     // Accessing nested values
 *     int subValue = static_cast<int>(std::get<jsonNumber>(std::get<ConfigJson::JsonMap>(config["sub"]).at("one")));
 *     std::cout << "Sub Value: " << subValue << std::endl;
 *
 *     return 0;
 * }
 * ```
 *
 * Example 4: Handling Errors
 *
 * ```cpp
 * #include "config_json.hpp"
 * #include <iostream>
 *
 * int main() {
 *     ConfigJson config;
 *     if (!config.parseJson("nonexistent.json")) {
 *         std::cerr << "Failed to parse JSON file." << std::endl;
 *         return 1;
 *     }
 *
 *     // Continue with JSON operations if parsing is successful
 *     std::cout << "Name: " << std::get<std::string>(config["name"]) << std::endl;
 *
 *     return 0;
 * }
 * ```
 *
 */

#ifndef CONFIG_JSON_HPP
#define CONFIG_JSON_HPP

#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

class jsonNumber {
public:
  jsonNumber(int value) : intValue(value), isInt(true) {}
  jsonNumber(double value) : doubleValue(value), isInt(false) {}

  operator int() const {
    if (isInt) {
      return intValue;
    } else {
      return static_cast<int>(doubleValue);
    }
  }

  operator double() const {
    if (isInt) {
      return static_cast<double>(intValue);
    } else {
      return doubleValue;
    }
  }

  // Overload the stream insertion operator
  friend std::ostream &operator<<(std::ostream &os, const jsonNumber &num) {
    if (num.isInt) {
      os << num.intValue;
    } else {
      os << num.doubleValue;
    }
    return os;
  }

private:
  int intValue;
  double doubleValue;
  bool isInt;
};

class ConfigJson {
public:
  struct JsonValue; // Forward declaration
  using JsonMap = std::map<std::string, JsonValue>;

  struct JsonValue : std::variant<std::string, jsonNumber, JsonMap> {
    using variant::variant;

    JsonValue &operator[](const std::string &key) {
      if (!std::holds_alternative<JsonMap>(*this)) {
        *this = JsonMap{};
      }
      auto &map = std::get<JsonMap>(*this);
      return map[key];
    }
  };

  ConfigJson() {}
  ConfigJson(const std::string &filepath) { parseJson(filepath); }

  void parseJson(const std::string &filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
      std::cout << "Could not open file: " << filepath << std::endl;
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();

    size_t pos = 0;
    jsonData = parseObject(content, pos);
  }

  JsonMap parseObject(const std::string &content, size_t &pos) {
    JsonMap object;
    skipWhitespace(content, pos);

    if (pos >= content.size() || content[pos] != '{') {
      std::cout << "Expected '{' at start of JSON object." << std::endl;
      return object; // Return an empty object if the start of the object is not found
    }

    pos++;

    while (pos < content.size()) {
      skipWhitespace(content, pos);

      if (content[pos] == '}') {
        pos++;
        return object;
      }

      std::string key = parseString(content, pos);
      skipWhitespace(content, pos);

      if (pos >= content.size() || content[pos] != ':') {
        std::cout << "Expected ':' after key in JSON object." << std::endl;
        return object; // Return the current state of the object if a colon is not found
      }

      pos++;
      JsonValue value = parseValue(content, pos);
      object[key] = value;
      skipWhitespace(content, pos);

      if (pos < content.size() && content[pos] == ',') {
        pos++;
      }
    }

    std::cout << "Expected '}' at end of JSON object." << std::endl;
    return object; // Return the current state of the object if the end is not found
  }

  void save(const std::string &filepath) {
    std::ofstream file(filepath);
    if (!file.is_open()) {
      std::cout << "Could not open file for writing: " << filepath << std::endl;
    }
    file << toString(jsonData);
    file.close();
  }

  JsonValue &operator[](const std::string &key) { return jsonData[key]; }

  const JsonValue &operator[](const std::string &key) const {
    static const JsonValue defaultValue = std::string("");
    auto it = jsonData.find(key);
    if (it != jsonData.end()) {
      return it->second;
    }
    return defaultValue;
  }

  template <typename T> T get(const std::string &key, const T defaultReturn = T{0}) {
    auto it = jsonData.find(key);
    if (it != jsonData.end()) {
      return std::get<T>(it->second);
    }
    return defaultReturn; // Default value if key is missing
  }

  template <typename T> T get_in(const std::string &key1, const std::string &key2, const T defaultReturn = T{0}) {
    auto it1 = jsonData.find(key1);
    if (it1 != jsonData.end()) {
      auto &innerMap = std::get<JsonMap>(it1->second);
      auto it2 = innerMap.find(key2);
      if (it2 != innerMap.end()) {
        return std::get<T>(it2->second);
      } 
    }
    return defaultReturn; // Default value if either key is missing
  }

private:
  JsonMap jsonData;

  std::string parseString(const std::string &content, size_t &pos) {
    skipWhitespace(content, pos);
    if (pos >= content.size() || content[pos] != '"') {
      std::cout << "Expected '\"' at start of string." << std::endl;
    }
    pos++;

    std::string str;
    while (pos < content.size() && content[pos] != '"') {
      str += content[pos++];
    }
    if (pos >= content.size()) {
      std::cout << "Expected '\"' at end of string." << std::endl;
    }
    pos++;
    return str;
  }

  JsonValue parseValue(const std::string &content, size_t &pos) {
    skipWhitespace(content, pos);
    if (pos >= content.size()) {
      std::cout << "Unexpected end of input while parsing value." << std::endl;
    }

    if (content[pos] == '"') {
      return parseString(content, pos);
    } else if (content[pos] == '{') {
      return parseObject(content, pos);
    } else {
      size_t endPos = pos;
      bool isDouble = false;
      while (endPos < content.size() &&
             (isdigit(content[endPos]) || content[endPos] == '.' || content[endPos] == '-')) {
        if (content[endPos] == '.') {
          isDouble = true;
        }
        endPos++;
      }

      std::string valueStr = content.substr(pos, endPos - pos);
      pos = endPos;

      if (isDouble) {
        return jsonNumber(std::stod(valueStr));
      } else {
        return jsonNumber(std::stoi(valueStr));
      }
    }
  }

  void skipWhitespace(const std::string &content, size_t &pos) {
    while (pos < content.size() && isspace(content[pos])) {
      pos++;
    }
  }

  std::string toString(const JsonMap &jsonMap, int indent = 1) {
    std::string jsonString;
    std::string indentStr(indent * 2, ' '); // 2 spaces per indent level

    jsonString += "{\n";
    for (auto it = jsonMap.begin(); it != jsonMap.end(); ++it) {
      if (it != jsonMap.begin()) {
        jsonString += ",\n";
      }
      jsonString += indentStr + "\"" + it->first + "\": " + valueToString(it->second, indent);
    }
    jsonString += "\n" + std::string((indent - 1) * 2, ' ') + "}";
    return jsonString;
  }

  std::string valueToString(const JsonValue &value, int indent = 1) {
    if (std::holds_alternative<std::string>(value)) {
      return "\"" + std::get<std::string>(value) + "\"";
    } else if (std::holds_alternative<jsonNumber>(value)) {
      jsonNumber num = std::get<jsonNumber>(value);
      if (static_cast<double>(static_cast<int>(num)) == static_cast<double>(num)) {
        return std::to_string(static_cast<int>(num));
      } else {
        return std::to_string(static_cast<double>(num));
      }
    } else if (std::holds_alternative<JsonMap>(value)) {
      return toString(std::get<JsonMap>(value), indent + 1);
    }
    return "";
  }
};

#endif // CONFIG_JSON_HPP

#if 0

#include <iostream>

int main() {
  try {
    // ConfigJson config("example.json");
    ConfigJson config;
    config.parseJson("example.json");
    /*
{
  "name": "toto tata",
  "time": 1.2,
  "age": 30,
  "height": 5.9,
  "sub": {
    "one": 1,
    "two": "test"
  }
}
    */

    std::cout << "int not found with default value = " << config.get<jsonNumber>("err", 12) << std::endl;
    std::string tt = config.get<std::string>("fake", "not found");
    std::cout << "string not found with default value = '" << tt << "'" << std::endl;
    std::cout << "error = " << config.get_in<jsonNumber>("sub", "err", 134) << std::endl;
    std::cout << "error = " << config.get_in<jsonNumber>("df", "err", 3.14159) << std::endl;

    int timeInt = std::get<jsonNumber>(config["time"]);
    double timeDouble = std::get<jsonNumber>(config["time"]);
    std::cout << "int = " << timeInt << ", double = " << timeDouble << std::endl;

    std::cout << "time = " << std::get<jsonNumber>(config["time"]) << std::endl;
    std::cout << "time = " << config.get<jsonNumber>("time") << std::endl;
    std::cout << "name = " << std::get<std::string>(config["name"]) << std::endl;
    std::cout << "age = " << std::get<jsonNumber>(config["age"]) << std::endl;
    std::cout << "height = " << std::get<jsonNumber>(config["height"]) << std::endl;
    std::cout << "sub one = "
              << std::get<jsonNumber>(std::get<std::map<std::string, ConfigJson::JsonValue>>(config["sub"])["one"])
              << std::endl;
    std::cout << "sub one = " << config.get_in<jsonNumber>("sub", "one") << std::endl;
    // std::cout << "sub one = " << config.get<std::string>(config["sub"]["two"]) << std::endl;

    // Save the configuration to a new file
    config.save("output.json");
    std::cout << "Configuration saved to output.json" << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}

#endif

#if 0

#include <iostream>

int main() {
  try {
    // ConfigJson config("example.json");
    ConfigJson config;
    config["param1"] = 0;
    config["param2"] = 0.2345;
    config["param3"] = "coucou";
    config["param3"] = "coucou mec";
    config["vector"]["x"] = 12.2;
    config["vector"]["y"] = 2.2;
    config["vector"]["z"] = -5.28;

    // Save the configuration to a new file
    config.save("output.json");
    std::cout << "Configuration saved to output.json" << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}

#endif

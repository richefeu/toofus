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

#ifndef STRING_MANIP_HPP
#define STRING_MANIP_HPP

#include <iostream>
#include <sstream>
#include <string>

#include <algorithm>
#include <cctype>
#include <codecvt>
#include <locale>
#include <unordered_set>
#include <vector>

inline void inplaceTrim(std::string &line) {
  // trim from start (in place)
  line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) { return !std::isspace(ch); }));
  // trim from end (in place)
  line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(),
             line.end());
}

inline std::string readQuotedString(std::istream &is) {
  char currentChar;
  std::string result;
  bool inQuotes  = false;
  char quoteChar = '\0';

  // Skip initial whitespace
  while (is.get(currentChar) && isspace(currentChar)) {}

  // Check if the input starts with a quote
  if (currentChar == '"' || currentChar == '\'') {
    inQuotes  = true;
    quoteChar = currentChar;
  } else {
    // If not quoted, add the first character back to the stream
    is.putback(currentChar);
  }

  if (inQuotes) {
    // Read until the closing quote or end of line
    std::getline(is, result, quoteChar);

    // Check if the line ended before finding the closing quote
    if (is.fail() && is.eof()) { std::cerr << "Warning: End-of-line reached before closing quote." << std::endl; }
  } else {
    // Read until the next space or end of line
    std::getline(is, result, ' ');
  }

  return result;
}

// Function to remove accents from a string
inline std::string removeAccents(const std::string &str) {
  std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> converter;
  std::u32string utf32 = converter.from_bytes(str);

  for (char32_t &c : utf32) {
    if ((c & 0xFF80) == 0x00C0) { // Check if it's a character with an accent
      c &= 0xFFDF;                // Remove the accent
    }
  }

  return converter.to_bytes(utf32);
}

// Function to normalize the string
inline std::string normalizeString(const std::string &str) {
  std::string normalized = removeAccents(str);
  std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::tolower);

  // Replace hyphens with spaces
  std::replace(normalized.begin(), normalized.end(), '-', ' ');

  // Remove extra spaces
  normalized.erase(
      std::unique(normalized.begin(), normalized.end(), [](char a, char b) { return a == ' ' && b == ' '; }),
      normalized.end());

  // Trim leading and trailing spaces
  normalized.erase(0, normalized.find_first_not_of(' '));
  normalized.erase(normalized.find_last_not_of(' ') + 1);

  return normalized;
}

// Function to split a string into words
inline std::unordered_set<std::string> tokenizeString(const std::string &str) {
  std::unordered_set<std::string> words;
  std::istringstream iss(str);
  std::string word;
  while (iss >> word) { words.insert(word); }
  return words;
}

// Function to check if a word starts with a given prefix
inline bool startsWith(const std::string &word, const std::string &prefix) {
  if (prefix.size() > word.size()) { return false; }
  return std::equal(prefix.begin(), prefix.end(), word.begin());
}

inline bool match(const std::string &s1, const std::string &s2, int &score) {
  std::string normalized1 = normalizeString(s1);
  std::string normalized2 = normalizeString(s2);

  std::unordered_set<std::string> words1 = tokenizeString(normalized1);
  std::unordered_set<std::string> words2 = tokenizeString(normalized2);

  // Count common words
  int commonWords = 0;
  for (const auto &word : words1) {
    if (words2.find(word) != words2.end()) { commonWords++; }
  }

  // Calculate similarity score
  int totalWords = std::max(words1.size(), words2.size());
  score          = (totalWords == 0) ? 0 : (commonWords * 100) / totalWords;

  // Consider it a match if at least half of the words match
  return commonWords >= totalWords / 2;
}

#endif // STRING_MANIP_HPP

#if 0
int main() {
    std::string input1 = "     Hello World     ";
    std::string input2 = "\"Quoted Text\" with the rest ";
    std::string input3 = "'Single Quoted'";
    std::string input4 = "\"Incomplete Quote \n ";

    std::istringstream iss1(input1);
    std::istringstream iss2(input2);
    std::istringstream iss3(input3);
    std::istringstream iss4(input4);

    std::string result1 = readQuotedString(iss1);
    std::string result2 = readQuotedString(iss2);
    std::string result3 = readQuotedString(iss3);
    std::string result4 = readQuotedString(iss4);

    std::cout << "Result 1: " << result1 << std::endl;
    std::cout << "Result 2: " << result2 << std::endl;
    std::cout << "Result 3: " << result3 << std::endl;
    std::cout << "Result 4: " << result4 << std::endl;

    return 0;
}
#endif

#if 0
int main() {
  std::string str1 = "Café--au-lait coucou";
  std::string str2 = "Lait au caféman";

  int score;
  
  bool result = match(str1, str2, score);
  std::cout << "Match: " << (result ? "True" : "False") << std::endl;
  std::cout << "Score: " << score << std::endl;
  
  result = match(str1, str1, score);
  std::cout << "Match: " << (result ? "True" : "False") << std::endl;
  std::cout << "Score: " << score << std::endl;

  return 0;
}
#endif

#include "doctest.h"
#include "fileTool.hpp"

TEST_CASE("fileTool") {
  REQUIRE(fileTool::fileExists("fileTool_test.cpp"));
  std::string str("/complicated/path/to/find/the/file/file.name.ext");
  REQUIRE(fileTool::GetFileExt(str) == "ext");
  REQUIRE(fileTool::GetFilePath(str) == "/complicated/path/to/find/the/file");
  REQUIRE(fileTool::GetFileName(str) == "file.name");
}

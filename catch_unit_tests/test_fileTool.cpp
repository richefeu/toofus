#include "catch.hpp"
#include "fileTool.hpp"

TEST_CASE("fileTool") {
  REQUIRE(fileTool::fileExists("test_fileTool.cpp"));
  std::string str("/complicated/path/to/find/the/file/file.name.ext");
  REQUIRE(fileTool::GetFileExt(str) == "ext");
  REQUIRE(fileTool::GetFilePath(str) == "/complicated/path/to/find/the/file");
  REQUIRE(fileTool::GetFileName(str) == "file.name");
}

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Werror -I./doctest -I..

OBJDIR = obj
SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))

TEST_EXEC = test_app

all: $(TEST_EXEC)

$(TEST_EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

run: $(TEST_EXEC)
	./$(TEST_EXEC) -s

clean:
	rm -rf $(OBJDIR) $(TEST_EXEC)

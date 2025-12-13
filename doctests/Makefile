CXX = g++
CXXFLAGS = -std=c++17 -Wall -Werror -I./doctest -I..

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

TEST_EXEC = test_app

all: $(TEST_EXEC)

$(TEST_EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TEST_EXEC)

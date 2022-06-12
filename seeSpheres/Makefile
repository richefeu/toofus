UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = clang++
  CXXFLAGS = -O3 -Wall -std=c++11
  GLFLAGS = `pkg-config --cflags --libs glut` -framework OpenGL 
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++0x
  GLFLAGS = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

all: seeSpheres

clean:
	rm -f *.o seeSpheres

seeSpheres: seeSpheres.cpp
	$(CXX) $(CXXFLAGS) -Wno-deprecated $(GLFLAGS) -o $@ $<

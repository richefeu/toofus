UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = clang++
  CXXFLAGS = -O3 -Wall -std=c++11
  
	#GLFLAGS = `pkg-config --cflags --libs glut` -framework OpenGL 
	
	# on apple, compile freeglut from source with this configuration:
	# mkdir build
	# cd build
	# cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -DFREEGLUT_COCOA=ON -DFREEGLUT_BUILD_DEMOS=OFF -DCMAKE_BUILD_TYPE=Release
  # make
	# sudo make install
	
  GLLINK = /usr/local/lib/libglut.a -framework Cocoa -framework CoreVideo -framework OpenGL -framework IOKit
	GLFLAGS = -I/usr/local/include
	
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++0x
  GLFLAGS = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

all: seeSpheres

clean:
	rm -f *.o seeSpheres

seeSpheres: seeSpheres.cpp
	$(CXX) $(CXXFLAGS) -Wno-deprecated $(GLFLAGS) $(GLLINK) -o $@ $<
	#$(CXX) $(CXXFLAGS) -c $< -o seeSpheres.o $(GLFLAGS)
	#$(CXX) -o $@ seeSpheres.o $(GLLINK)

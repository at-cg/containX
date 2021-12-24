EXE= containX
LIBS= -lz
#CPPFLAGS= -std=c++11 -g
#CPPFLAGS= -std=c++11 -O3
CPPFLAGS= -std=c++11 -O3 -fopenmp
#CPPFLAGS= -std=c++11 -O3 -DVERBOSE

all:
	$(CXX) $(CPPFLAGS) main.cpp -o $(EXE) $(LIBS) 

clean:
	rm -fr $(EXE)


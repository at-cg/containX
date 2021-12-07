EXE= containX
LIBS= -lz
CPPFLAGS= -std=c++11 -O3

all:
	$(CXX) $(CPPFLAGS) main.cpp -o $(EXE) $(LIBS) 

clean:
	rm -fr $(EXE)


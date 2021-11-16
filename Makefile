EXE=containX
LIBS=   -lz

all:
	$(CXX) $(CPPFLAGS) main.cpp -o $(EXE) $(LIBS) 

clean:
	rm -fr $(EXE)


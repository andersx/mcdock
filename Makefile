CXX = g++

OBDIR = /home/andersx/opt/openbabel2.4

INCLUDE = -I$(OBDIR)/include/openbabel-2.0/
LIBS = -L$(OBDIR)/lib

all: mcdock

mcdock: src/main_ob.cpp
	$(CXX) $(INCLUDE) $(LIBS) -std=c++11 -O3 -march=native src/main_ob.cpp -o mcdock -lopenbabel

conformers: src/conformers.cpp
	$(CXX) $(INCLUDE) $(LIBS) -std=c++11 -O3 -march=native src/conformers.cpp -o conformers -lopenbabel

clean:
	rm -f main



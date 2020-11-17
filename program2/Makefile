CXX = g++
CXXFLAGS = -Wall -pedantic -O3 -std=c++17
LFLAGS =

OBJS = main.o

all: program

program: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f *.o

.PHONY: all clean
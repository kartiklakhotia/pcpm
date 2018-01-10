CC      = g++
CPPFLAGS= -O3 -c -std=c++11 -fopenmp -w -mavx 
LDFLAGS = -fopenmp -m64 -lpthread 
SOURCES = pr.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) pr

pr : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o pr dump*


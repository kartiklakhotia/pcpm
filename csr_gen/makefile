#HUGE_EDGE=1
#HUGE_VERTEX=1

ifdef HUGE_EDGE
INTE = -DHUGE_EDGE
endif

ifdef HUGE_VERTEX
INTV = -DHUGE_VERTEX
endif

CC      = g++
CPPFLAGS= -O3 -c -std=c++11 -w $(INTE) $(INTV)
LDFLAGS = -static -O3 -m64
SOURCES = main.cpp graph.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) a.out

a.out : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o a.out


CC      = g++
PCM_PATH = /home/klakhoti/graphChallenge/pcm
CPPFLAGS= -O3 -c -std=c++11 -fopenmp -mavx -w 
LDFLAGS = -fopenmp -m64 -lpthread 
SOURCES = pr.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) pr

pr : $(OBJECTS)  

	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o pr dump*


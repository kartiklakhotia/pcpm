CC      = g++
PCM_PATH = /home/klakhoti/graphChallenge/pcm
CPPFLAGS= -O3 -c -std=c++11 -w
LDFLAGS = -m64
SOURCES = relabel.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) relabel

relabel : $(OBJECTS)  

	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o relabel dump*


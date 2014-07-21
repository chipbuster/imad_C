CC=g++
LD=g++

OBJS=imad.o GdalFileIO.o geo_utils.o imad_utils.o ImageStats.o imad_bigfun.o
LIBS=-lgdal -lm
INCLUDEDIRS= -I/usr/include/eigen3

CCFLAGS= -O3 -g -p -c -Wall
LDFLAGS=-g -p

all: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) $(LIBS) -o imad

imad.o: imad.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) imad.cpp -o imad.o

imad_utils.o: imad_utils.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) imad_utils.cpp -o imad_utils.o

geo_utils.o: geo_utils.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) geo_utils.cpp -o geo_utils.o

imad_bigfun.o: imad_bigfun.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) imad_bigfun.cpp -o imad_bigfun.o

GdalFileIO.o: GdalFileIO.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) GdalFileIO.cpp -o GdalFileIO.o

ImageStats.o: ImageStats.cpp
	$(CC) $(CCFLAGS) $(INCLUDEDIRS) ImageStats.cpp -o ImageStats.o

debug: imad
	gdb imad

test: all
	./imad


clean:
	rm *.o imad

CC=g++
LD=g++
CCFLAGS= -O3 -g -p -c -Wall -I/usr/include/eigen3
LDFLAGS=-g -p -lgdal -lm


all: imad.o GdalFileIO.o imad_utils.o ImageStats.o geo_utils.o
	$(CC) imad.o GdalFileIO.o imad_utils.o ImageStats.o $(LDFLAGS) -o imad

imad.o: imad.cpp
	$(CC) $(CCFLAGS) imad.cpp -o imad.o

imad_utils.o: imad_utils.cpp
	$(CC) $(CCFLAGS) imad_utils.cpp -o imad_utils.o

geo_utils.o:
	$(CC) $(CCFLAGS) geo_utils.cpp -o geo_utils.o

GdalFileIO.o: GdalFileIO.cpp
	$(CC) $(CCFLAGS) GdalFileIO.cpp -o GdalFileIO.o

ImageStats.o: ImageStats.cpp
	$(CC) $(CCFLAGS) ImageStats.cpp -o ImageStats.o

debug: imad
	gdb imad

profile:


clean:
	rm *.o imad

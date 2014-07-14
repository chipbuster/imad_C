CC=g++
LD=g++
CCFLAGS=-c -Wall -I/usr/include/eigen3
LDFLAGS=-lgdal -lm


all: imad.o GdalFileIO.o imad_utils.o ImageStats.o
	$(CC) imad.o GdalFileIO.o imad_utils.o ImageStats.o $(LDFLAGS) -o imad

imad.o: imad.cpp
	$(CC) $(CCFLAGS) imad.cpp -o imad.o

imad_utils.o: imad_utils.cpp
	$(CC) $(CCFLAGS) imad_utils.cpp -o imad_utils.o

GdalFileIO.o: GdalFileIO.cpp
	$(CC) $(CCFLAGS) GdalFileIO.cpp -o GdalFileIO.o

ImageStats.o: ImageStats.cpp
	$(CC) $(CCFLAGS) ImageStats.cpp -o ImageStats.o

clean:
	rm *.o imad

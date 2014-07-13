CC=g++
LD=g++
CCFLAGS=-c -Wall -I/usr/include/eigen3
LDFLAGS=-lgdal


all: imad.o GdalFileIO.o
	$(CC) imad.o GdalFileIO.o $(LDFLAGS) -o imad

imad.o: imad.cpp
	$(CC) $(CCFLAGS) imad.cpp -o imad.o

GdalFileIO.o: GdalFileIO.cpp
	$(CC) $(CCFLAGS) GdalFileIO.cpp -o GdalFileIO.o

clean:
	rm *.o imad

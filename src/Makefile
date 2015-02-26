CC=g++
CFLAGS=-Wall

all: main

main: Distance.o main.cc
	$(CC) $(CFLAGS) Distance.o main.cc -o main

Distance.o: Distance.h Distance.cc
	$(CC) $(CFLAGS) -c Distance.cc

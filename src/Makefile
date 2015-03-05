CC=g++
CFLAGS=-g -c -Wall

all: main

main: Distance.o main.o
	$(CC) Distance.o main.o -o main

main.o: Distance.o main.cc
	$(CC) $(CFLAGS) main.cc

Distance.o: Distance.h Distance.cc
	$(CC) $(CFLAGS) Distance.cc

run:
	./main

clean:
	rm main ./*.o

.PHONY: all run clean

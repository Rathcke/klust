CC=g++
CFLAGS=-g -c -Wall -std=c++11 -O3 -march=native #-pg
LDFLAGS= #-pg

all: main

main: Distance.o main.o Cluster.o IO.o
	$(CC) $(LDFLAGS) Cluster.o Distance.o IO.o main.o -o main

main.o: Distance.o main.cc
	$(CC) $(CFLAGS) main.cc

Distance.o: Distance.h Distance.cc
	$(CC) $(CFLAGS) Distance.cc

Cluster.o: Cluster.h Cluster.cc
	$(CC) $(CFLAGS) Cluster.cc

IO.o: IO.h IO.cc
	$(CC) $(CFLAGS) IO.cc

run:
	./main

clean:
	rm main ./*.o

.PHONY: all run clean

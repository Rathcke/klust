CC=g++
CFLAGS=-g -c -Wall -std=c++11 -O3 -march=native -pthread #-pg
#LDFLAGS=-pthread #-pg
LDFLAGS=-pthread -static # -static-libgcc -static-libstdc++ -static-glibc

all: main

main: Distance.o main.o Cluster.o IO.o Seq.o
	$(CC) $(LDFLAGS) Cluster.o Distance.o IO.o Seq.o main.o -o main

main.o: Distance.o main.cc
	$(CC) $(CFLAGS) main.cc

Distance.o: Distance.h Distance.cc
	$(CC) $(CFLAGS) Distance.cc

Cluster.o: Cluster.h Cluster.cc
	$(CC) $(CFLAGS) Cluster.cc

IO.o: IO.h IO.cc
	$(CC) $(CFLAGS) IO.cc

Seq.o: Seq.h Seq.cc
	$(CC) $(CFLAGS) Seq.cc

run:
	./main

clean:
	rm main ./*.o

.PHONY: all run clean

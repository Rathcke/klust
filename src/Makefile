CC=g++
CFLAGS=-g -c -Wall -std=c++11 -O3 -march=native -pthread
LDFLAGS=-pthread

all: klust

klust: Distance.o klust.o Cluster.o IO.o Seq.o Utils.o
	$(CC) $(LDFLAGS) Cluster.o Distance.o IO.o Seq.o Utils.o klust.o -o klust

klust.o: Distance.o klust.cc
	$(CC) $(CFLAGS) klust.cc

Distance.o: Distance.h Distance.cc
	$(CC) $(CFLAGS) Distance.cc

Cluster.o: Cluster.h Cluster.cc
	$(CC) $(CFLAGS) Cluster.cc

IO.o: IO.h IO.cc
	$(CC) $(CFLAGS) IO.cc

Seq.o: Seq.h Seq.cc
	$(CC) $(CFLAGS) Seq.cc

Utils.o: Utils.h Utils.cc
	$(CC) $(CFLAGS) Utils.cc

run:
	./klust

clean:
	rm klust ./*.o

.PHONY: all run clean

CC=g++
MSG=$(shell date)
STD=c++11
THREAD=NOTHREADS
libmd.o: libmd.cc libmd.h libmd/*
	$(CC) -Wall -Wextra -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -D $(THREAD) -pthread -fopenmp -std=$(STD) -O3 -c libmd.cc

all:
	make
	make -C ./tests

clean:
	rm libmd.o

cleanall:
	make clean
	make -C ./tests clean

forceclean:
	git clean -f

CC=g++
MSG=$(shell date)
libmd.o: libmd.cc libmd.h libmd/*
	$(CC) -Wall -Wextra -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -std=c++11 -O3 -c libmd.cc

all:
	make
	make -C ./test

clean:
	rm libmd.o

cleanall:
	make clean
	make -C ./test clean

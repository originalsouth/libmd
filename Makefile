CC=g++
MSG=$(shell date)
STD=c++11
libmd.o: libmd.cc libmd.h libmd/*
	$(CC) -Wall -Wextra -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -std=$(STD) -O3 -c libmd.cc

all:
	make
	make -C ./tests

clean:
	rm libmd.o

cleanall:
	make clean
	make -C ./tests clean

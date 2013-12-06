CC=g++
MSG=$(shell date)
STD=c++11
THREAD=NOTHREADS
curved-orbit: curved-orbit.cc ../../libmd.cc ../../libmd.h ../../libmd/*
	$(CC) -Wall -Wextra -Wno-varargs -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -D $(THREAD) -pthread -fopenmp -std=$(STD) -O3 -o curved-orbit curved-orbit.cc

clean:
	rm curved-orbit

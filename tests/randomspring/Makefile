CC=g++
MSG=$(shell date)
STD=c++11
THREAD=NOTHREADS
randomspring: springio.cc randomspring.cc ../../libmd.cc ../../libmd.h ../../libmd/*
	$(CC) -Wall -Wextra -Wno-varargs -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -D $(THREAD) -pthread -fopenmp -std=$(STD) -O3 -o randomspring randomspring.cc

clean:
	rm randomspring

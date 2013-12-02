CC=g++
MSG=$(shell date)
STD=c++11
THREAD=NOTHREADS
md: md.cc ../../libmd.cc ../../libmd.h ../../libmd/*
	$(CC) -Wall -Wextra -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -D $(THREAD) -pthread -fopenmp -std=$(STD) -O3 -o md md.cc

clean:
	rm md
	rm *.png

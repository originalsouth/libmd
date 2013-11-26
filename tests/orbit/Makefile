CC=g++
MSG=$(shell date)
md: md.cc ../../libmd.cc ../../libmd.h ../../libmd/*
	$(CC) -Wall -Wextra -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -std=c++11 -O3 -o md md.cc

clean:
	rm md


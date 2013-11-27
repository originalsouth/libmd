CC=g++
MSG=$(shell date)
STD=c++11
md: md.cc ../../libmd.cc ../../libmd.h ../../libmd/*
	$(CC) -Wall -Wextra -D CMSG='"$(MSG)"' -D CC='"$(CC)"' -std=$(STD) -O3 -o md md.cc

clean:
	rm md
	rm *.png

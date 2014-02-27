CC=g++
MSG=$(shell date)
VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
STD=c++11
THREAD=NOTHREADS
CCWFLAGS=-Wall -Wextra 
CCDFLAGS=-D CMSG='"$(MSG)"' -D CC='"$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL)"' -D $(THREAD) -D VER='"$(VER)"' -D BRANCH='"$(BRANCH)"'
CCOPTLEVEL=3
freelyjointedchain: freelyjointedchain.cc freelyjointedchain2.cc ../../libmd.cc ../../libmd.h ../../libmd-src/*
ifeq ($(THREAD),OPENMP)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -fopenmp $(CCDFLAGS) -o freelyjointedchain freelyjointedchain.cc
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -fopenmp $(CCDFLAGS) -o freelyjointedchain2 freelyjointedchain2.cc
else ifeq ($(THREAD),THREADS)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -pthread $(CCDFLAGS) -o freelyjointedchain freelyjointedchain.cc
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -pthread $(CCDFLAGS) -o freelyjointedchain2 freelyjointedchain2.cc
else
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) $(CCDFLAGS) -o freelyjointedchain freelyjointedchain.cc
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) $(CCDFLAGS) -o freelyjointedchain2 freelyjointedchain2.cc
endif
clean:
	rm freelyjointedchain
	rm freelyjointedchain2

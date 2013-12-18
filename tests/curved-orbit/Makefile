CC=g++
MSG=$(shell date)
VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
STD=c++11
THREAD=NOTHREADS
CCWFLAGS=-Wall -Wextra -Wno-varargs
CCDFLAGS=-D CMSG='"$(MSG)"' -D CC='"$(CC)"' -D $(THREAD) -D VER='"$(VER)"' -D BRANCH='"$(BRANCH)"'
CCOPTLEVEL=3
curved-orbit: curved-orbit.cc ../../libmd.cc ../../libmd.h ../../libmd/*
ifeq ($(THREAD),OPENMP)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -fopenmp $(CCDFLAGS) -o curved-orbit curved-orbit.cc
else ifeq ($(THREAD),THREADS)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -pthread $(CCDFLAGS) -o curved-orbit curved-orbit.cc
else
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) $(CCDFLAGS) -o curved-orbit curved-orbit.cc
endif

clean:
	rm curved-orbit

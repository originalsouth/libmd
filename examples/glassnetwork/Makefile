CC=g++
MSG=$(shell date)
VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
STD=c++11
THREAD=NOTHREADS
CCWFLAGS=-Wall -Wextra 
CCDFLAGS=-D CMSG='"$(MSG)"' -D CC='"$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL)"' -D $(THREAD) -D VER='"$(VER)"' -D BRANCH='"$(BRANCH)"'
CCOPTLEVEL=3
glass-shear: glass-shear.cc ../../libmd.cc ../../libmd.h ../../libmd-src/* ../../tools/springio/springio.cc ../../tools/springio/PointSystem2d/*
ifeq ($(THREAD),OPENMP)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -fopenmp $(CCDFLAGS) -o glass-shear glass-shear.cc
else ifeq ($(THREAD),THREADS)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -pthread $(CCDFLAGS) -o glass-shear glass-shear.cc
else
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) $(CCDFLAGS) -o glass-shear glass-shear.cc
endif

clean:
	rm glass-shear

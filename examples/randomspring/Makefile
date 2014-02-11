CC=g++
MSG=$(shell date)
VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
STD=c++11
THREAD=NOTHREADS
CCWFLAGS=-Wall -Wextra -Wno-varargs
CCDFLAGS=-D CMSG='"$(MSG)"' -D CC='"$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL)"' -D $(THREAD) -D VER='"$(VER)"' -D BRANCH='"$(BRANCH)"'
CCOPTLEVEL=3
randomspring: randomspring.cc ../../libmd.cc ../../libmd.h ../../libmd-src/* ../../tools/springio/springio.cc ../../tools/springio/PointSystem2d/*
ifeq ($(THREAD),OPENMP)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -fopenmp $(CCDFLAGS) -o randomspring randomspring.cc
else ifeq ($(THREAD),THREADS)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -pthread $(CCDFLAGS) -o randomspring randomspring.cc
else
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) $(CCDFLAGS) -o randomspring randomspring.cc
endif

clean:
	rm randomspring

CC=g++
MSG=$(shell date)
VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
STD=c++11
THREAD=NOTHREADS
CCWFLAGS=-Wall -Wextra
CCDFLAGS=-D CMSG='"$(MSG)"' -D CC='"$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL)"' -D $(THREAD) -D VER='"$(VER)"' -D BRANCH='"$(BRANCH)"'
CCOPTLEVEL=3
libmd.o: libmd.cc libmd.h libmd-src/*
ifeq ($(THREAD),OPENMP)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -fopenmp $(CCDFLAGS) -c libmd.cc
else ifeq ($(THREAD),THREADS)
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) -pthread $(CCDFLAGS) -c libmd.cc
else
	$(CC) $(CCWFLAGS) -std=$(STD) -O$(CCOPTLEVEL) $(CCDFLAGS) -c libmd.cc
endif
	

all:
	make
	make tests
	make examples
	make doc

.PHONY: tests
tests:
	make -C ./tests

.PHONY: examples
examples:
	make -C ./examples

.PHONY: doc 
doc:
	make -C ./doc

clean: 
	rm libmd.o

.PHONY: cleantests
cleantests:
	make -C ./tests clean

.PHONY: cleanexamples
cleanexamples: 
	make -C ./examples clean

.PHONY: cleandoc
cleandoc: 
	make -C ./doc clean

cleanall:
	make clean
	make cleantests
	make cleanexamples
	make cleandoc

forceclean:
	git clean -f

include Makeheader

libmd.o: libmd.cc libmd.h libmd-src/*
	$(CC) $(CCFLAGS) -c libmd.cc

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

include Makeheader

.PHONY: libmd.o

#PLUGDIR=$(shell gcc -print-search-dirs | awk -F\: '/install/{print $$2}')
#libmd.so: libmd.cc libmd-src/*.cc libmd-src/md/*.cc libmd.h
#	mkdir -p obj
#	for i in `find . -name "*.cc" | grep libmd-src`; do $(CC) $(CCFLAGS) -ffat-lto-objects -fPIC -o obj/`basename $$i | sed "s/\.cc/\.o/g"` -c $$i; done
#	gcc-ar -rcs libmd.a `find obj -name "*.o" -printf "%p "` --plugin$(PLUGDIR)liblto_plugin.so
#	$(CC) -flto -shared -fPIC libmd.a -o libmd.so 

libmd.so: libmd.cc libmd-src/*.cc libmd-src/md/*.cc libmd.h
	$(CC) -c $(CCFLAGS) -fPIC -o libmd.o libmd.cc
	$(CC) $(CCFLAGS) -shared -fPIC -o libmd.so libmd.o

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

#clean: 
#	rm -rf obj
#	rm libmd.a libmd.so

clean: 
	rm -rf obj
	rm libmd.o libmd.so

.PHONY: cleantests
cleantests:
	make -C ./tests clean

.PHONY: cleanexamples
cleanexamples: 
	make -C ./examples clean

.PHONY: cleandoc
cleandoc: 
	make -C ./doc cleanall

cleanall:
	make clean
	make cleantests
	make cleanexamples
	make cleandoc

forceclean:
	git clean -f

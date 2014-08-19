VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

.PHONY: all 
all:
	make intro
	make doxygen

.PHONY: doxygen
doxygen:
	(sed -u -e "/PROJECT_NUMBER         =/ s/= .*/= $(BRANCH)-0.$(VER)/" Doxyfile) | doxygen -

.PHONY: intro
intro:
	make -C ./libmd-intro

cleanall:
	make cleandoxygen
	make cleanintro

cleandoxygen:
	rm -rf doxygen

.PHONY: cleanintro
cleanintro:
	make -C ./libmd-intro clean

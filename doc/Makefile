VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
doxygenVERtest=$(awk -F'.' '{print ($2<8?1:0)}' <(doxygen --version)) #Remove line if Roland ever updates the system

.PHONY: all 
all:
	make intro
	make doxygen

doxygen:
ifeq ($(doxygenVERtest),1) #Remove line if Roland ever updates the system
	(sed -u -e "/PROJECT_NUMBER         =/ s/= .*/= $(BRANCH)-0.$(VER)/" Doxyfile) | doxygen -std=c++0x - #Remove if Roland ever updates the system
else #Remove line if Roland ever updates the system
	(sed -u -e "/PROJECT_NUMBER         =/ s/= .*/= $(BRANCH)-0.$(VER)/" Doxyfile) | doxygen -
endif #Remove line if Roland ever updates the system

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

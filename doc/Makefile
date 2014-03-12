VER=$(shell git rev-list HEAD --count)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
doxygenVERtest=$(awk -F'.' '{print ($2<8?1:0)}' <(doxygen --version)) #Remove line if Roland ever updates the system
Doxyfile-doc: ../*
	sed -e "/PROJECT_NUMBER         =/ s/= .*/= $(BRANCH)-0.$(VER)/" Doxyfile > Doxyfile-doc
ifeq ($(doxygenVERtest),1) #Remove line if Roland ever updates the system
	doxygen -std=c++0x Doxyfile-doc #Remove if Roland ever updates the system
else #Remove line if Roland ever updates the system
	doxygen Doxyfile-doc
endif #Remove line if Roland ever updates the system

clean:
	rm -rf ./html ./latex ./man ./xml Doxyfile-doc

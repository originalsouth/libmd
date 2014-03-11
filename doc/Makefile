VER=$(shell git rev-list HEAD --count)
doxygenVERtest=$(awk -F'.' '{print ($2<8?1:0)}' <(doxygen --version))
Doxyfile-doc: ../*
	sed -e "/PROJECT_NUMBER         =/ s/= .*/= 0.$(shell git rev-list HEAD --count)/" Doxyfile > Doxyfile-doc
ifeq ($(doxygenVERtest),1)
	doxygen -std=c++0x Doxyfile-doc #Remove if Roland ever updates system
else
	doxygen Doxyfile-doc
endif

clean:
	rm -rf ./html ./latex ./man ./xml Doxyfile-doc

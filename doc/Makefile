VER=$(shell git rev-list HEAD --count)
doc: ../*
	sed -i -e "/PROJECT_NUMBER         =/ s/= .*/= 0.$(shell git rev-list HEAD --count)/" Doxyfile
	touch .doc
	doxygen Doxyfile
	
clean:
	rm -rf ./html ./latex ./man ./xml .doc

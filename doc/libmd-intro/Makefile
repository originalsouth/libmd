.PHONY: libmd-intro.pdf all clean

all: libmd-intro.pdf

%.tex: %.raw
	./raw2tex $< > $@

%.tex: %.dat
	./dat2tex $< > $@

libmd-intro.pdf: libmd-intro.tex meta.xmp
	latexmk -pdf -quiet -interaction=nonstopmode -use-make libmd-intro.tex

clean:
	latexmk -CA
	rm meta.xmpi

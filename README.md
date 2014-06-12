libmd -- A molecular dynamics library optimized for soft matter 
=======

### Some stuff libmd does:
  * It aims for nice 'black-box' like (but highly hackable) interface making md as easy (or hard) as you want
  * It runs in any integer dimension (but gets somewhat slower)
  * It can take care of Newton's equation on Riemannian geometry

### Required/recommended software:
  * gcc (4.8.x or higher), alternatively clang (3.3 or higher), or a C++11 complaint compiler
  * git
  * GNU Make
  * awk/sed 
  * Doxygen
  * LaTeX/LaTeX-Mk
  * Mogrify (ImageMagick) for libmd-projects using BaX (see tools/BaX)

### Using libmd in your project:
  * The easiest way to start a new libmd project is via the pre-built libmd-empty-project hosted [here](https://bitbucket.org/softmatter/empty-libmd-project)
  * Another way is to develop your project in the example tree
  * Finally you can compile libmd as a library/shared-library and install it in your system

### Learning libmd:
  * See doc/libmd-intro
  * There a dozen of examples and tests in examples/ and test/
  * Bother and author (see AUTHORS file)

### Developing libmd: 
  * For now there is doc/doxygen to help you (don't forget to make it)

### Developers
The developers of libmd are listed in the AUTHORS file

Last updated: Thu Jun 12 11:59:46 CEST 2014

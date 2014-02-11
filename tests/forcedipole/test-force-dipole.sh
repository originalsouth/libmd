#!/bin/bash

# requires gnuplot
NEWCC='/data/misc/simulshare/compilers/gcc-4.8.2/bin/g++'
make CC=$NEWCC forcedipole

echo "Running forcedipole..."
./forcedipole > __foo
echo "...done."

grep pos __foo > __bar

gnuplot ./test-force-dipole-graph.gp
rm __foo __bar
echo "See test-force-dipole.eps for numerics vs. theory"

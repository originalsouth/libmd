set key center right
set xlabel 'Frame'
set ylabel 'x-coordinate'

set term postscript enhanced color
set output 'test-force-dipole.eps'
plot '__bar' u ($2) ti 'Point 5', '' u ($4) ti 'Point 6', -(.45+.05 *cos(sqrt(2)*.1*x)) notitle, (.45+.05*cos(sqrt(2)*.1*x)) notitle
unset output


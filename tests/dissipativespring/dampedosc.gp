#!/usr/bin/gnuplot

set term postscript enhanced color
set size 0.7, 0.7
set samples 500
set xlabel 't'
set ylabel 'x(t)'
set datafile separator " "

# initial conditions
v0 = 0.1
x0 = 0 
w0 = 1

# dissipative coefficient = 0.1 underdamped
xi = 0.1
wd = w0*sqrt(1-xi*xi)
a = x0
b = (v0 + w0*x0*xi)/wd

set output 'underdamp.eps'
plot '< grep pos ./output_c0.1' u 2:($9-1) ti 'Sim', exp(-xi*w0*x)*(a*cos(wd*x)+b*sin(wd*x)) ti 'Theory'


# dissipative coefficient = 1 critically damped
a = x0
b = v0 + w0*x0
set output 'critdamp.eps'
plot '< grep pos ./output_c1.0' u 2:($9-1) ti 'Sim', (a+b*x)*exp(-w0*x) ti 'Theory'

# dissipative coefficient = 10 overdamped
xi = 10.
gp = (-2.*xi*w0+sqrt(4.*xi*xi*w0*w0-4.*w0*w0))/2.
gm = (-2.*xi*w0-sqrt(4.*xi*xi*w0*w0-4.*w0*w0))/2.
a = x0 - (gm*x0-v0)/(gm-gp)
b = (gm*x0-v0)/(gm-gp)

set output 'overdamp.eps'
plot '< grep pos ./output_c10.0' u 2:($9-1) ti 'Sim', a*exp(gm*x)+b*exp(gp*x) ti 'Theory'


#!/usr/bin/gnuplot
set terminal wxt persist
set datafile separator ";"
file="<./autodiff"
plot file u 1:2 w l title "f1",\
     file u 1:3 w l title "autodiff f1",\
     file u 1:4 w l title "finitediff f1",\
     file u 1:5 w l title "f2",\
     file u 1:6 w l title "autodiff f2",\
     file u 1:7 w l title "finitediff f2"

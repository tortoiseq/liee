#------------------------------------------------------------------------------
load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL
set output template.serial.".png"

maIN=data; @EXTREMA
dk = (k1 - k0) / (Nk - 1)

set xrange [k0_plt:k1_plt]
set yrange [glob_min:glob_max]
set xlabel "k in 1/nm"
set ylabel "probability density"

plot data using (k0 + $0*dk):($1) with lines title "detection rate"

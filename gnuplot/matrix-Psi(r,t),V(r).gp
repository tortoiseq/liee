# example plot of the wf's Psi^2-matrix ontop of an overview of the min/max potential
#------------------------------------------------------------------------------

load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL

set output template.serial.".png"
set multiplot

maIN=data; @EXTREMA

set lmargin at screen r0_scr
set rmargin at screen r1_scr
set bmargin at screen 0.90
set tmargin at screen 0.4
load '/usr/local/src/gnuplot/sub_matrix-Psi(r,t).gp'

set lmargin at screen r0_scr
set rmargin at screen r1_scr
set bmargin at screen 0.4
set tmargin at screen 0.1
load '/usr/local/src/gnuplot/sub_V(r)_@ttt.gp'

load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL

set output template.serial.".png"
set multiplot

@SCAN_DATA

# calculate potential data
system( sprintf( "liee_tools.py potential %f %f %d %f %f %d 'potential.dat'", r0, r1, Nr, t0, t1, Nt) )
system( "grep '##' 'potential.dat' | awk '{print $2}' > 'potential.dat.param'" )
load 'potential.dat.param'

set lmargin at screen r0_scr
set rmargin at screen r1_scr
set bmargin at screen 0.95
set tmargin at screen 0.4
load '/usr/local/src/gnuplot/sub_matrix-Psi(r,t).gp'

set lmargin at screen r0_scr
set rmargin at screen r1_scr
set bmargin at screen 0.4
set tmargin at screen 0.05
load '/usr/local/src/gnuplot/sub_V(r)_@ttt.gp'

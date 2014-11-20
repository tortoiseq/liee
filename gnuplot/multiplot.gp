# plots the Psi(r,t) as matrix-image and alligns with it V(r) and V(t), j(t) to 
# frame a comprehansive picture of the experiment
#------------------------------------------------------------------------------

load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL

set output template.serial.".png"
set multiplot

set lmargin at screen r0_scr
set rmargin at screen r1_scr
set bmargin at screen 0.90
set tmargin at screen 0.4
data=data_Psi
load '/usr/local/src/gnuplot/sub_matrix-Psi(r,t).gp'

set lmargin at screen r0_scr
set rmargin at screen r1_scr
set bmargin at screen 0.4
set tmargin at screen 0.1
load '/usr/local/src/gnuplot/sub_V(r)_@ttt.gp'

set lmargin at screen r1_scr
set rmargin at screen 0.9
set bmargin at screen 0.9
set tmargin at screen 0.4
data=data_Jt
load '/usr/local/src/gnuplot/sub_V(t),J(t).gp'

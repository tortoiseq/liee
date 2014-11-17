load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL
set output template.serial.".png"

load '/usr/local/src/gnuplot/sub_V(r)_@t.gp'

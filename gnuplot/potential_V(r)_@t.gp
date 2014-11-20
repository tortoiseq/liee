# simple example for just the potential V(r) at a given time.
# mainly for testing if the program routines work
#------------------------------------------------------------------------------


load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL
set output template.serial.".png"

load '/usr/local/src/gnuplot/sub_V(r)_@t.gp'

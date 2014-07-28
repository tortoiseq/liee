load 'set_vars.gp'
set term png enhanced size 1366, 768
set output template.serial.".png"

Nr_plt = 1000
Nt = 1

system( sprintf( "liee_tools.py potential %f %f %d %f %f %d 'potential.dat'", r0_plt, r1_plt, Nr_plt, t_plt, t_plt, Nt) )

set xrange [r0_plt:r1_plt]
set xlabel "r in nm"
set xtics
set ytics
set yrange [V0_plt:V1_plt]
set ylabel "potential in eV"

plot 'potential.dat' using 1:2 with lines title "V(r)"

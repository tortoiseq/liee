Nt = 1
system( sprintf( "liee_tools.py potential %f %f %d %f %f %d 'potential.dat'", r0_plt, r1_plt, Nr_plt, t_plt, t_plt, Nt) )
system( "grep '##' 'potential.dat' | awk '{print $2}' > 'potential.dat.param'" )
load 'potential.dat.param'

set xrange [r0_plt:r1_plt]
set xlabel "r in nm"
set xtics
set ytics
set yrange [V_min - 0.03*(V_max - V_min) : V_max + 0.03*(V_max - V_min)]
set ylabel "potential in eV"

plot 'potential.dat' using 1:2 with lines title "V(r)"

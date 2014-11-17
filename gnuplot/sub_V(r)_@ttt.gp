set xrange [r0_plt:r1_plt]
set xlabel "r in nm"
set xtics
set ytics
set yrange [V_min - 0.03*(V_max - V_min) : V_max + 0.03*(V_max - V_min)]
set ylabel "potential in eV"
unset key

mini=int(V_min_t+2)
maxi=int(V_max_t+2)

plot 'potential.dat' using 1:2 with lines title "V_0",\
     'potential.dat' using 1:mini with lines title "V_{min}",\
     'potential.dat' using 1:maxi with lines title "V_{max}"

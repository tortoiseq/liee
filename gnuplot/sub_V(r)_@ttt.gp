# for plotting three V(r) curves: at start and +/- maximum Laser Field
#------------------------------------------------------------------------------

# calculate potential data
maIN="V(r,t).dat"; ma_r0=r0_plt; ma_r1=r1_plt; ma_Nr=Nr_plt; ma_t0=t0; ma_t1=t1; ma_Nt=Nt; @POTENTIAL
dr = (r1 - r0) / (Nr - 1)

set xrange [r0_plt:r1_plt]
set xlabel "r in nm"
set xtics
set ytics
set yrange [Vmin - 0.03*(Vmax - Vmin) : Vmax + 0.03*(Vmax - Vmin)]
set ylabel "potential in eV"
set ytics nomirror
unset key
set style line 1 linetype  1 linewidth 2 linecolor rgb "blue"
set style line 2 linetype  1 linewidth 2 linecolor rgb "black"
set style line 3 linetype  1 linewidth 2 linecolor rgb "dark-magenta"

mini=int(Vmin_i+2)
maxi=int(Vmax_i+2)

plot 'V(r,t).dat' using (r0+$0*dr):2 with lines ls 2 title "V_0",\
     'V(r,t).dat' using (r0+$0*dr):mini with lines ls 1 title "V_{min}",\
     'V(r,t).dat' using (r0+$0*dr):maxi with lines ls 3 title "V_{max}"

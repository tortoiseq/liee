# example 101 of sub-plot
#------------------------------------------------------------------------------

maIN="V(r).dat"; ma_r0=r0_plt; ma_r1=r1_plt; ma_Nr=Nr_plt; ma_t0=t_plt; ma_t1=t_plt; ma_Nt=1; @POTENTIAL
dr = (r1 - r0) / (Nr - 1)

set xrange [r0_plt:r1_plt]
set xlabel "r in nm"
set xtics
set ytics
set yrange [Vmin - 0.03*(Vmax - Vmin) : Vmax + 0.03*(Vmax - Vmin)]
set ylabel "potential in eV"

plot 'V(r).dat' using (r0+$0*dr):1 with lines title "V(r)"

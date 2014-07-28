# select closest data record for plot_t of this frame
t_idx = 1 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt - 1) )

set xrange [r0_plt:r1_plt]
set yrange [glob_min:glob_max]
set xlabel "r in nm"
set xtics
set ylabel "probability density"
set ytics nomirror
set y2tics
set y2range [-0.2:11]
set y2label "potential in eV"

plot 'wf(r,t)_.dat' using (r0 + $0*dr):t_idx with lines title sprintf( "{/Symbol Y}(r, %.3g fs)", t_plt ),\
     'potential.dat' using 1:2 with lines axes x1y2 title "V(r)"


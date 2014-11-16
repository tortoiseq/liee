# select closest data record for plot_t of this frame
t_idx = 1 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt - 1) )
pt_idx = 2 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt_plt - 1) )

set xrange [r0_plt:r1_plt]
set yrange [glob_min:glob_max]
set xlabel "r in nm"
set xtics
set ylabel "probability density"
set ytics nomirror
set y2tics
set y2range [-200:200]
set y2label "potential in eV"

plot data.".t" using (r0 + $0*dr):t_idx with lines title sprintf( "{/Symbol Y}(r, %.3f fs)", t_plt ),\
     'potential.dat' using 1:pt_idx with lines axes x1y2 title "V(r)"


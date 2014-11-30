# incomplete
#------------------------------------------------------------------------------


# select closest data record for plot_t of this frame
t_idx = floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt - 1) )

set xrange [r0_plt : r1_plt]
set yrange [k0 : k1]
set cbrange [glob_max/10.0 : glob_max]
#unset colorbox
unset y2tics
unset y2label
set ylabel "k in 1/nm"
set xlabel "r in nm"

plot data using (r0 + $2*dr):(k0 + $1*dk):3 matrix index t_idx with image notitle

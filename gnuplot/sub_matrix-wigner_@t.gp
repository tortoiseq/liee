# incomplete
#------------------------------------------------------------------------------


# select closest data record for plot_t of this frame
#t_idx = 1 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt - 1) )
dr = (r1 - r0) / (Nr - 1)

#print( sprintf( "t_idx=%d   t_plt=%f   t0=%f    t1=%f    Nt=%d", t_idx, t_plt, t0, t1, Nt) )

#set xrange [r0_plt:r1_plt]
#set yrange [1e3*hbar/me * k0 : 1e3*hbar/me * k1]
#set cbrange [0:glob_max]
unset colorbox
unset xtics
unset xlabel
unset y2tics
unset y2label
#set ylabel "v in nm/fs"
set ylabel "~ v in nm/fs"
#plot data using (r0 + $2*dr):(1e3*hbar/me*(k0 + $1*dk)):3 matrix index t_idx with image
plot data using 2:1:3 matrix index frame with image

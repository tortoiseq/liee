# Template for plotting V(t) rotated 90 degrees to align with a matrix plot besides
#
# required variables: 	key_x, key_y, Vmin, Vmax, tVmin, tVmax, dt, Nt, t_charge
# required data-files:	data - current over time
#------------------------------------------------------------------------------


maIN=data; @PARAMETER
dt = (t1 - t0) / (Nt - 1)
maIN="V(t).dat"; ma_r0=r_detect; ma_r1=r_detect; ma_Nr=1; ma_t0=t0; ma_t1=t1; ma_Nt=Nt; @POTENTIAL
maIN="V(t).dat"; @TRANSPOSE

maIN=data; @EXTREMA
set style line 1 linetype  1 linewidth 1 linecolor rgb "blue" 
set style line 2 linetype  1 linewidth 1 linecolor rgb "red" 

#set key at screen key_x, key_y
unset ytics  
unset ylabel 
set y2tics 
set y2label 'Zeit in fs' 
set y2range [0: Nt*dt]

set format x "%.1tx10^{%2T}"
set xtics 
set xlabel 'Emissionsrate in 1/fs' #offset 0,-2.5 
set xrange[glob_min:glob_max]
set xtics rotate 
set x2tics rotate 
set x2label 'Potential in eV' #offset 0,2 
set x2range[Vmin-1e9:Vmax+1e9]

#set obj 1 circle at second Vmin*CONVeV, Vmin_t*dt*CONVfs  size screen 0.008 front fillcolor rgb "blue" fillstyle solid 1.0 
#set obj 2 circle at second Vmax*CONVeV, Vmax_t*dt*CONVfs  size screen 0.008 front fillcolor rgb "dark-magenta" fillstyle solid 1.0 
#set obj 3 circle at second tc_V*CONVeV, t_charge  size screen 0.008 front fillcolor rgb "black" fillstyle solid 1.0 

plot "V(t).dat.t" using 1:($0*dt) axis x2y2 title 'Pulse-Potential' with lines ls 1,\
     data using 1:($0*dt) axis x1y2 title 'Tunnel Strom' with lines ls 2


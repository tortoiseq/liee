# incomplete
#------------------------------------------------------------------------------

if(plot_frame==-1 && !skip_scan){
    system( sprintf( "grep '##' 'P(r,k).dat' | awk '{print $2}' > 'P(r,k).param'" ) )
}
load 'P(r,k).param'
dr = (r1 - r0) / (Nr - 1)
dk = (k1 - k0) / (Nk - 1)

# scan through data to find global extrema
if(plot_frame==-1 && !skip_scan){
  print "scan through P(r,k)-data to find global extrema"
  system( sprintf( "python2.7 /home/quark/workspace/LIEEpy/src/liee/tools/liee_tools.py stats 'P(r,k).dat' 'P(r,k).param'") )
  load 'P(r,k).param'
}
else{
  # select closest data record for plot_t of this frame
  t_idx = floor( 0.5 + (plot_t - t0) / (t1 - t0) * (Nt - 1) )

  set xrange [plot_r0:plot_r1]
  set yrange [1e3*hbar/me*k0:1e3*hbar/me*k1]
  set cbrange [0:glob_max]
  unset colorbox
  unset xtics
  unset xlabel
  unset y2tics
  unset y2label
  set ylabel "v in nm/fs"
  plot 'P(r,k).dat' using (r0 + $2*dr):(1e3*hbar/me*(k0 + $1*dk)):3 matrix index t_idx with image
  
}


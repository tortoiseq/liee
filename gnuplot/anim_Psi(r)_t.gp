load 'set_vars.gp'
set term png enhanced size 1366, 768

# extract some parameters from the datafile
system( sprintf( "grep '##' 'wf(r,t).dat' | awk '{print $2}' > 'wf(r,t).param'" ) )
load 'wf(r,t).param'
dr = (r1 - r0) / (Nr - 1)

# scan through data to find global extrema
system( sprintf( "liee_tools.py stats 'wf(r,t).dat' 'wf(r,t).param'" ) )
load 'wf(r,t).param'

# calculate potential data
system( sprintf( "liee_tools.py potential %f %f %d %f %f %d 'potential.dat'", r0, r1, Nr, t0, t1, Nt) )

# transpose data so that gnuplot can operate on columns instead of rows
system( sprintf( "liee_tools.py transpose 'wf(r,t).dat' 'wf(r,t)_.dat'" ) )
system( sprintf( "grep '#' 'wf(r,t).dat' >> 'wf(r,t)_.dat'" ) )


dt_plt = (t1_plt - t0_plt) / (Nt_plt - 1)


do for [frame=0 : Nt_plt-1]{
  t_plt = t0_plt + frame * dt_plt
  print "frame ".frame
  set output sprintf( "%s.%d.%03d.png", template, serial, frame )

  set lmargin at screen r0_scr
  set rmargin at screen r1_scr
  set bmargin at screen 0.05
  set tmargin at screen 0.95
  load '/usr/local/src/gnuplot/sub_snapshot_Psi(r).gp'
}


system( sprintf( "mencoder 'mf://./%s'.%d.???.png -mf type=png:w=1366:h=768:fps=5 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o '%s'.%d.avi", template, serial, template, serial ) )

# clean up
system( sprintf( "rm '%s'.%d.???.png", template, serial ) )
system( sprintf( "rm 'wf(r,t)_.dat'" ) )

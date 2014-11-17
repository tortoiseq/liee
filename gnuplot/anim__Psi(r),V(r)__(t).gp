# Creats an animation of the 1d-probability density together with the time-varying potential.
#
# @calls( sub_snapshot_Psi(r).gp ) for each frame
# @param( data ) -> data-file with wave-function samples, e.g. from Obs_Snapshot_WF
# @param( r0, r1, Nr, t0, t1, Nt ) -> space-time sampling grid of the recorded data is auto-extracted from data-file
# @param( r0_plt, r1_plt, t0_plt, t1_plt ) -> portion of the data-range to plot
# @param( Nt_plt ) -> total frames (or time resolution)
# @param( width_plt, height_plt, fps_plt ) -> resolution of plot, speed of animation
#
# TODO auto-detect data-range of the potential
#------------------------------------------------------------------------------------------------------------------------------

load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL
@SCAN_DATA

# calculate potential data
system( sprintf( "liee_tools.py potential %f %f %d %f %f %d 'potential.dat'", r0, r1, Nr, t0, t1, Nt_plt) )
system( "grep '##' 'potential.dat' | awk '{print $2}' > 'potential.dat.param'" )
load 'potential.dat.param'

# transpose data so that gnuplot can operate on columns instead of rows
system( sprintf( "liee_tools.py transpose '%s' '%s.t'", data, data ) )
system( sprintf( "grep '#' '%s' >> '%s.t'", data, data ) )

dt_plt = ( t1_plt - t0_plt ) / ( Nt_plt - 1 )

do for [frame=0 : Nt_plt-1]{
  t_plt = t0_plt + frame * dt_plt
  print "frame ".frame
  set output sprintf( "%d.%05d.png", serial, frame )
  load '/usr/local/src/gnuplot/sub_snapshot_Psi(r),V(r)_@t.gp'
}


system( sprintf( "mencoder 'mf://./%d.?????.png' -mf type=png:w=%d:h=%d:fps=%d -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o '%s.%d.avi'", serial, width_plt, height_plt, fps_plt, template, serial ) )

# clean up
system( sprintf( "rm %d.?????.png", serial ) )
system( sprintf( "rm '%s.t'", data ) )

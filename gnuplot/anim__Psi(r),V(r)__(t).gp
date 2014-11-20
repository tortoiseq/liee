# Creats an animation of the 1d-probability density together with the time-varying potential.
#
# @calls( sub_snapshot_Psi(r),V(r)_@t.gp ) for each frame
# @param( data ) -> data-file with wave-function samples, e.g. from Obs_Snapshot_WF
# @param( r0, r1, Nr, t0, t1, Nt ) -> space-time sampling grid of the recorded data is auto-extracted from data-file
# @param( r0_plt, r1_plt, t0_plt, t1_plt ) -> portion of the data-range to plot
# @param( Nt_plt ) -> total frames (or time resolution)
# @param( width_plt, height_plt, fps_plt ) -> resolution of plot, speed of animation
#
#------------------------------------------------------------------------------

load 'set_vars.gp'
load '/usr/local/src/gnuplot/liee_macros.gp'
@TERMINAL

maIN=data; @PARAMETER
maIN='V(r,t).dat'; ma_r0=r0; ma_r1=r1; ma_Nr=Nr; ma_t0=t0; ma_t1=t1; ma_Nt=Nt_plt; @POTENTIAL
data_V='V(r,t).dat'
dt_plt = ( t1_plt - t0_plt ) / ( Nt_plt - 1 )

maIN=data; @EXTREMA; @TRANSPOSE
data=data.".t"

do for [frame=0 : Nt_plt-1]{
  t_plt = t0_plt + frame * dt_plt
  set output sprintf( "%d.%05d.png", serial, frame )
  load '/usr/local/src/gnuplot/sub_snapshot_Psi(r),V(r)_@t.gp'
}

system( sprintf( "mencoder 'mf://./%d.?????.png' -mf type=png:w=%d:h=%d:fps=%d -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o '%s.%d.avi'", serial, width_plt, height_plt, fps_plt, template, serial ) )
system( sprintf( "rm %d.?????.png", serial ) )  # clean up

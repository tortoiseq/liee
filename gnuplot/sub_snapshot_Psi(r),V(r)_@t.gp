# Plots an instance in time of the wave-function together with the current potential
# (a sub-plot is called from a super-plot which sets up terminal and margins and provides parameters)
#
# @precondition( file: potential.dat ) grided potential data according to plot-range and resolution
# @param ( t_plt ) -> time requested to plot. the actual plot is from the nearest available time in the data (no interpolation)
# @param( r0, r1, Nr, t0, t1, Nt ) -> space-time sampling grid of the recorded data is auto-extracted from data-file
# @param( r0_plt, r1_plt ) -> portion of the data-range to plot
# @param( Nt_plt ) -> total frames (or time resolution)
# @param( glob_min, glob_max ) -> global bounds of the wave-function data (for scaling)
#------------------------------------------------------------------------------------------------------------------------------

# select closest data record for plot_t of this frame
t_idx = 1 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt - 1) )
pt_idx = 2 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt_plt - 1) )
dr = (r1 - r0) / (Nr - 1)


set xrange [r0_plt:r1_plt]
set yrange [glob_min:glob_max]
set xlabel "r in nm"
set xtics
set ylabel "probability density"
set ytics nomirror
set y2tics
set y2range [V_min - 0.03*(V_max - V_min) : V_max + 0.03*(V_max - V_min)]
set y2label "potential in eV"

plot data.".t" using (r0 + $0*dr):t_idx with lines title sprintf( "{/Symbol Y}(r, %.3f fs)", t_plt ),\
     'potential.dat' using 1:pt_idx with lines axes x1y2 title "V(r)"


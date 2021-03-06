# Plots an instance in time of the wave-function together with the current potential
# (a sub-plot is called from a super-plot which sets up terminal and margins and provides parameters)
#
# @precondition( file with filename in variable: data ) grided wave-function data
# @precondition( file with filename in variable: data_V ) grided potential data according to plot-range and resolution
# @param ( t_plt ) -> time requested to plot. the actual plot is from the nearest available time in the data (no interpolation)
# @param( r0, r1, Nr, t0, t1, Nt ) -> space-time sampling grid of the recorded data (should be auto-extracted from data-file)
# @param( r0_plt, r1_plt ) -> portion of the data-range to plot
# @param( Nt_plt ) -> total frames (or time resolution)
# @param( glob_min, glob_max ) -> global bounds of the wave-function data (for scaling)
#------------------------------------------------------------------------------

# select closest data record for plot_t of this frame
t_idx = 1 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt - 1) )
pt_idx = 1 + floor( 0.5 + (t_plt - t0) / (t1 - t0) * (Nt_plt - 1) )
dr = (r1 - r0) / (Nr - 1)

set xrange [r0_plt:r1_plt]
set yrange [glob_min:glob_max]
set xlabel "r in nm"
set xtics
set ylabel "probability density"
set ytics nomirror
set y2tics
set y2range [Vmin - 0.03*(Vmax - Vmin) - 1e6 : Vmax + 0.03*(Vmax - Vmin) + 1e6]
set y2label "potential in eV"

plot data using (r0 + $0*dr):t_idx with lines title sprintf( "{/Symbol Y}(r, %.3f fs)", t_plt ),\
     data_V using (r0 + $0*dr):pt_idx with lines axes x1y2 title "V(r)"


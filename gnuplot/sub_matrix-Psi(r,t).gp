# (sub-plot to be called from another gnuplot file, for setting up the stage)
# plots two-dimensional data {e.g. wave function squares Psi(r,t)^2 or its Fourier-
# transformed Psi(k,t)^2 } as matrix-image-plot with a logarithmic color-scale in grey. 
# the scale is truncated to [log10(Psi_huge) : log10(Psi_tiny)]
#------------------------------------------------------------------------------

maIN=data; @EXTREMA

dr = (r1 - r0) / (Nr - 1)
dt = (t1 - t0) / (Nt - 1)

#unset xtics 
#unset ytics 
set xrange [r0_plt:r1_plt]
set yrange [t0_plt:t1_plt]

set cbrange [log10(Psi_tiny):log10(Psi_huge)] 
#set colorbox horizontal user origin (d1+cbar_hborder),(sz-d4+cbar_vborder) size (sz-m-d1-2*cbar_hborder), (d4-2*cbar_vborder) front 
unset colorbox
unset key
#set label "log_{10}( {/Symbol-Oblique Y}^2 nm )" at screen 0.055*sz, screen 0.95*z front
set pm3d map 
set palette defined( log10(Psi_tiny) 1 1 1, (log10(Psi_huge)-0.03) 0 0 0, log10(Psi_huge) 1 0.8 0.8 ) 

plot data matrix using (r0+$1*dr):(t0+$2*dt):($3<Psi_tiny?log10(Psi_tiny):($3>Psi_huge?log10(Psi_huge):log10($3)) ) with image




# for Psi(k,t): conversion to nm/fs
#set ylabel "v in nm/fs"
#set yrange [1e3*hbar/me*k0:1e3*hbar/me*k1]
#plot 'P(r,k).dat' using (r0 + $2*dr):(1e3*hbar/me*(k0 + $1*dk)):3 matrix index t_idx with image

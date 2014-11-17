unset xtics 
unset ytics 
set xrange [0:Nr]
set yrange [0:Nt]

#set cbrange [log10(TINY):log10(HUGE)] 
#set colorbox horizontal user origin (d1+cbar_hborder),(sz-d4+cbar_vborder) size (sz-m-d1-2*cbar_hborder), (d4-2*cbar_vborder) front 
unset colorbox
#set label "log_{10}( {/Symbol-Oblique Y}^2 nm )" at screen 0.055*sz, screen 0.95*z front
set pm3d map 
#set palette defined( log10(TINY) 1 1 1, (log10(HUGE)-0.03) 0 0 0, log10(HUGE) 1 0.8 0.8 ) 

#plot data matrix using 1:2:($3<TINY?log10(TINY):($3>HUGE?log10(HUGE):log10($3)) ) with image
plot data matrix using 1:2:3 with image

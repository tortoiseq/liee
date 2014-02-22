'''
Created on 04-Feb-2012

Tool for plotting the simulation snapshot files taken by liee::Obs_Snapshot_WF.

@author: quark
'''
import numpy
import sys
import subprocess
import getopt
import math
import ctypes

CONV_au_m = 5.2917720859e-11
CONV_au_nm = 5.2917720859e-2
CONV_au_eV = 27.211
CONV_au_s = 2.418884326505e-17
CONV_au_fs = 2.418884326505e-2
CONV_au_V_over_m = 5.1421e11
CONV_au_V = CONV_au_V_over_m * CONV_au_m

def main():
    # parse arguments, load datafile, find maxima of Psi, minima of V
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:pt:r:l", ["file=","potential","tiny=","huge=","logscale","movie","complex"])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    TINY = float('nan')
    HUGE = float('nan')
    filename = ""
    do_pot = False
    logscale = False
    do_square = True
    do_movie = False
    t_range = False 
    r_range = False
    
    for o, a in opts:
        if o in ("-f", "--file"):
            filename = a
        elif o in ("-p", "--potential"):
            do_pot = True
        elif o in ("-t",):
            t_range = float( a )/CONV_au_fs
        elif o in ("-r",):
            r_range = float( a )/CONV_au_nm
        elif o in ("--tiny"):
            TINY = float( a )
        elif o in ("--huge"):
            HUGE = float( a )
        elif o in ("-l", "--logscale"):
            logscale = True
        elif o in ("--movie"):
            do_movie = True
        elif o in ("--complex"):
            do_square = False
        else:
            assert False, "unhandled option"

#####################################################################################
#                       Acquire Data                                                #
#####################################################################################

    j_tot = 0
    E0 = 0
    r_0 = 0
    depth = 1
    width = 1
    t_charge = 0

    fobj = open("summary.txt", "r")
    for line in fobj: 
        if ( r_range == False ) and ( line.find( "r_range\t" ) > -1 ):
            i = line.find( "r_range\t" )
            r_range = float( line[( i + len("r_range\t") ):] ) /CONV_au_nm 
        elif ( line.find( "r_start\t" ) > -1 ):
            i = line.find( "r_start\t" )   
            r_0 = float( line[( i + len("r_start\t") ):] ) /CONV_au_nm
        elif ( t_range == False ) and ( line.find( "t_range\t" ) > -1 ):
            i = line.find( "t_range\t" )   
            t_range = float( line[( i + len("t_range\t") ):] ) /CONV_au_fs
        elif ( line.find( "init_E0\t" ) > -1 ):
            i = line.find( "init_E0\t" )   
            E0 = float( line[( i + len("init_E0\t") ):] )
        elif ( line.find( "tunnel_ratio\t" ) > -1 ):
            i = line.find( "tunnel_ratio\t" )   
            j_tot = float( line[( i + len("tunnel_ratio\t") ):] )
        elif ( line.find( "depth\t" ) > -1 ):
            i = line.find( "depth\t" )   
            depth = float( line[( i + len("depth\t") ):] )
        elif ( line.find( "width\t" ) > -1 ):
            i = line.find( "width\t" )   
            width = float( line[( i + len("width\t") ):] )
        elif ( line.find( "t_charge\t" ) > -1 ):
            i = line.find( "t_charge\t" )   
            t_charge = float( line[( i + len("t_charge\t") ):] ) 
    fobj.close()
    
    if ( do_square == True ):
        psi = numpy.loadtxt( filename, numpy.float, '#', )
    else:
        num_col = len( open( filename ).readline().split( "\t" ) )
        psi = numpy.loadtxt( filename, numpy.float, comments='#', usecols = range( 0, num_col-1, 2 ) )
        psi_im = numpy.loadtxt( filename, numpy.float, comments='#', usecols = range( 1, num_col-1, 2 ) )

    Nt = psi.shape[0]
    Nr = psi[0].shape[0]
    print(Nt)
    print(Nr)
    psi_log = numpy.zeros( (Nt, Nr) )
    dr = r_range / ( Nr - 1 )
    dt = t_range / ( Nt - 1 )
    print( "t_charge = " + str(t_charge) )
    print( "t_range = " + str(t_range) )
    print( "Nt = " + str(Nt) )
    tc_i = int( t_charge / CONV_au_fs / dt ) + 1 
    print( "tc_i = " + str(tc_i) )
    
    # auto-scale margins
    min = 1e+66
    max = 1e-66
    for ti in range( Nt ):
        for ri in range( Nr ):
            absr = abs( psi[ti][ri] )
            min = absr if ( absr < min ) else min
            max = absr if ( absr > max ) else max
            if ( not do_square ):
                psi[ti][ri] = psi[ti][ri] / CONV_au_nm
                absi = abs( psi_im[ti][ri] )
                min = absi if ( absi < min ) else min
                max = absi if ( absi > max ) else max
    TINY = min if ( math.isnan( TINY ) ) else TINY
    HUGE = max if ( math.isnan( HUGE ) ) else HUGE
    print( "Bounds:  [ " + str( TINY ) + " : " + str( HUGE ) + " ]\n")
    
    pot = None
    if do_pot:
        lib = ctypes.cdll.LoadLibrary('/usr/local/lib/libliee.so')
        lib.init_potential()
        pot = numpy.zeros( (Nt, Nr) )
        c_Vr = ctypes.c_double( 0.0 );
        c_Vi = ctypes.c_double( 0.0 );
        for ti in range( Nt ):
            for ri in range( Nr ):
                c_r = ctypes.c_double( r_0 + ri * dr );
                c_t = ctypes.c_double( ti * dt );
                lib.calc_potential( c_r, c_t, ctypes.byref(c_Vr), ctypes.byref(c_Vi) )
                pot[ti][ri] = c_Vr.value

    if ( not do_square ):
        outfile = open( "matrix.dat", "w" )
    for ti in range( Nt ):
        for ri in range( Nr ):
            if psi[ti][ri] < TINY:
                psi_log[ti][ri] = math.log10( TINY )
            else:
                psi_log[ti][ri] = math.log10( psi[ti][ri] )
            if ( not do_square ):
                #TODO logscale support
                realb = int( 255 * ( abs( psi[ti][ri] ) - TINY ) / (HUGE - TINY) )
                imagr = int( 255 * ( abs( psi_im[ti][ri] ) - TINY ) / (HUGE - TINY) )
                outfile.write( str(ri) + "\t" + str(ti) + "\t" + str( imagr ) + "\t" + str(0) + "\t" + str( realb ) + "\n" )

    if ( do_square ):
        if ( logscale ):
            numpy.savetxt( "matrix.dat", psi_log )
        else:
            numpy.savetxt( "matrix.dat", psi )
    else:
        outfile.close()
        
    r_split = 50.0
    i_split = int( ( r_split / CONV_au_nm - r_0 ) / dr )
    print( str(i_split) + " of "+ str(Nr) )
    tVmin = 0
    tVmax = 0
    if ( do_pot ):
        Vmin = pot[0][Nr-1]
        Vmax = pot[0][Nr-1]
        V_of_t = numpy.zeros( Nt )
        V_of_r = numpy.zeros( Nr )
        for i in range( Nt ):
            V_of_t[i] = pot[i][i_split]
            if V_of_t[i] < Vmin:    
                Vmin = V_of_t[i]
                tVmin = i
            if V_of_t[i] > Vmax:    
                Vmax = V_of_t[i]
                tVmax = i
        numpy.savetxt( "pulse_shape.dat", V_of_t )
    
        tc_V = pot[tc_i][i_split]
        for i in range( Nr ):
            V_of_r[i] = pot[tc_i][i]
        numpy.savetxt( "V(t0,r).dat", V_of_r )

        Vmin_r = 0.0
        Vmax_r = -0.0
        Vmin_l = 0.0
        Vmax_l = -0.0
        for i in range( Nr ):
            V_of_r[i] = pot[tVmin][i]
            if V_of_r[i] < Vmin_l and ( (r_0 + i*dr)*CONV_au_nm < r_split ) and ( (r_0 + i*dr)*CONV_au_nm > -width/2 ):
                Vmin_l = V_of_r[i]
            if V_of_r[i] < Vmin_r:    
                Vmin_r = V_of_r[i]
            if ( V_of_r[i] > Vmax_r ) and ( r_0 + i*dr > 0.0 ):
                Vmax_r = V_of_r[i]
        numpy.savetxt( "V(r)_pull.dat", V_of_r )

        for i in range( Nr ):
            V_of_r[i] = pot[tVmax][i]
            if V_of_r[i] > Vmax_l and ( (r_0 + i*dr)*CONV_au_nm < r_split ) and ( (r_0 + i*dr)*CONV_au_nm > -width/2 ):
                Vmax_l = V_of_r[i]
            if V_of_r[i] < Vmin_r:
                Vmin_r = V_of_r[i]
            if ( V_of_r[i] > Vmax_r ) and ( r_0 + i*dr > 0.0 ):
                Vmax_r = V_of_r[i]
        numpy.savetxt( "V(r)_push.dat", V_of_r )
        
        dV = Vmax_r - Vmin_r
        Vmax_r += 0.1 * dV
        Vmin_r -= 0.1 * dV


#####################################################################################
#                   Potential slices                                                #
#####################################################################################
    if ( do_pot ):
        gp = "set term png size 800, 600\n"
        gp += "set ylabel  \"V(r)\"\n"
        #gp += "set yrange [" + str( Vmin ) + ":" + str( Vmax ) + "]\n"
        gp += "set xlabel 'r (in arbitrary units)'\n"
        gp += "#set key off\n"
        gp += "#t = sprintf(\"%5.3g\",time)\n"
        gp += "set output 'V_of_r.png'\n"
        gp += "plot \"-\" using (" + str(r_0) + " + $0 * " + str(dr) + "):1 title \"V_min\" with lines, '' using (" + str(r_0) + " + $0 * " + str(dr) + "):1 title \"V_max\" with lines \n"
        stream_in = gp
        for x in range( Nr ):
            stream_in += str( pot[tVmin][x] ) + "\n"
        stream_in += "e\n"
        for x in range( Nr ):
            stream_in += str( pot[tVmax][x] ) + "\n"
        stream_in += "e\n"
        proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
        print proc.communicate( stream_in )[0]
    
    
        gp = "set term png size 800, 600\n"
        gp += "set ylabel  \"V(t)\"\n"
        gp += "#set yrange [-1.4:0.1]\n"
        gp += "set xlabel \"t (in arbitrary units)\"\n"
        gp += "#set key off\n"
        gp += "#t = sprintf(\"%5.3g\",time)\n"
        gp += "set output 'V_of_t.png'\n"
        gp += "plot \"-\" using ($0 * " + str(dt) + "):1 title \"Puls-Potential\" with lines \n"
        stream_in = gp
        for i in range( Nt ):
            stream_in += str( V_of_t[i] ) + "\n"
        stream_in += "e\n"
        proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
        print proc.communicate( stream_in )[0]
    
#####################################################################################
#                   Matrix plot                                                     #
#####################################################################################    
    
    gnuplot0 = "set term png size 1200, 800\n"
    gnuplot0 += "set xrange [0:" + str( Nr - 1 ) + "]\n"
    gnuplot0 += "set yrange [0:" + str( Nt - 1 ) + "]\n"
    gnuplot0 += "set output 'matrix.png'\n"
    if ( do_square ):
        gnuplot0 += "set pm3d map\n"
        gnuplot0 += "plot \"matrix.dat\" matrix with image\n"
    else:
        gnuplot0 += "plot \"matrix.dat\" with rgbimage\n"

    proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
    print proc.communicate( gnuplot0 )[0]

#####################################################################################
#                     MULTIPLOT                                                     #
#####################################################################################    
    #gp =  "set terminal epslatex color solid 'Helvetica' 20 \n"
    gp =  "set terminal postscript enhanced eps color solid 'Helvetica' 20 \n"
    #gp =  "set terminal cairolatex eps color solid font 'Helvetica' fontscale 1 \n"
    gp += "set output 'multiplot.ps' \n"
    gp += "set key off \n"

    gp += "sz = 2.5 \n"
    gp += "m = 0.2 * sz \n"
    gp += "d1 = 0.038 * sz \n" 
    gp += "d2 = 0.055 * sz \n"  # matching margin for t-axis (bigger d2 -> downwards)
    gp += "d3 = 0.03 * sz \n"   # horizontaly moves the V(t)-plot
    gp += "d4 = 0.09 * sz \n"   # matching margin for j-axis of t-plot
    gp += "blowup = 1.15 \n"    # blows up V-plots to make better use of the canvas and close the gap between V(r) and matrix plot
    gp += "cbar_hborder = 0.25 \n"
    gp += "cbar_vborder = 0.05 \n"

    gp += "set size sz, sz \n"
    gp += "set origin 0, 0 \n"
    gp += "set multiplot \n"

    # r-axis
    if ( do_pot ):
#        # backup
#        gp += "set style line 1 linetype  1 linewidth 3 linecolor rgb \"blue\" \n"
#        gp += "set size sz, m * blowup \n"
#        gp += "set lmargin 8.2 \n"
#        gp += "set rmargin 33 \n"
#        gp += "set origin 0, 0 \n"
#        gp += "set ylabel 'Potential in eV' \n"
#        gp += "set xlabel 'Entfernung von der Oberflaeche in nm' \n"
#        gp += "set yrange [" + str( Vmin_r * CONV_au_eV ) + ":" + str( Vmax_r * CONV_au_eV ) + "]\n"
#        gp += "set xrange [" + str( r_0 * CONV_au_nm ) + ":" + str( (r_0 + Nr * dr) * CONV_au_nm ) + "]\n"
#        gp += "plot 'V(r)_pull.dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_min\" with lines ls 1," 
#        gp +=       "'V(r)_push.dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_max\" with lines ls 1\n"

        r_range = Nr * dr * CONV_au_nm
        sz = 158.2  # whats the size of the plot in character width?
        x_range = sz - 33.0 - 8.2
        x_split = 8.2 + ( r_split - r_0 * CONV_au_nm ) * x_range / r_range

        stretch = abs( Vmin_l )
        stretch *= CONV_au_eV * 1.07

        # Potential-Well
        gp += "set style line 1 linetype  1 linewidth 4 linecolor rgb \"blue\" \n"
        gp += "set style line 2 linetype  1 linewidth 4 linecolor rgb \"black\" \n"
        gp += "set style line 3 linetype  1 linewidth 4 linecolor rgb \"dark-magenta\" \n"
        gp += "set size sz, m * blowup \n"
        gp += "set lmargin 8.2 \n"
        gp += "set rmargin " + str( 33 + x_range - x_split ) + " \n"
        gp += "set bmargin 4.0 \n"
        gp += "set origin 0, 0 \n"
        gp += "set ylabel 'Potential in eV' \n"
        gp += "set ytics nomirror\n"
        gp += "unset xlabel \n"
        gp += "set yrange [" + str( -1.0 * stretch ) + ":" + str( 1.0 * stretch ) + "]\n"
        gp += "set xrange [" + str( r_0 * CONV_au_nm ) + ":" + str( r_split ) + "]\n"
        gp += "plot 'V(r)_pull.dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_min\" with lines ls 1,"
        gp +=      "'V(r)_push.dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_max\" with lines ls 3,"
        gp +=      "'V(t0,r).dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_0\" with lines ls 2\n"

        # Potential-Tail
        stretch = abs( Vmin_r )
        if abs( Vmax_r ) > stretch:
            stretch = abs( Vmax_r )
        stretch *= CONV_au_eV * 1.05
        
        gp += "set size sz, m * blowup \n"
        gp += "set lmargin " + str( x_split ) + " \n"
        gp += "set rmargin 33 \n"
        gp += "set bmargin 4.0 \n"
        gp += "set origin 0, 0 \n"
        gp += 'set format y "%.1tx10^{%2T}" \n'
        gp += "yRange=" + str( stretch ) + "\n"
        gp += "set ytics -0.75*yRange, 0.375*yRange, 0.75*yRange \n"
        gp += "set ytics mirror offset " + str( 117.0 - x_split ) + ", 0 \n"
        gp += "unset ylabel \n"
        gp += "set xlabel 'Entfernung von der Oberflaeche in nm' \n"
        gp += "set yrange [-yRange:yRange]\n"
        gp += "set xrange [" + str( r_split ) + ":" + str( (r_0 + Nr * dr) * CONV_au_nm ) + "]\n"
        gp += "pos(x) = x < " + str( 0.90 * Nr ) + " ? " + str( r_0 * CONV_au_nm ) + " + x * " + str( dr * CONV_au_nm ) + " : NaN \n"
        gp += "pot(x) = x * " + str( CONV_au_eV ) + " \n"
        gp += "plot 'V(r)_pull.dat' using (pos($0)):(pot($1)) title \"V_min\" with lines ls 1, " 
        gp +=        "'V(t0,r).dat' using (pos($0)):(pot($1)) title \"V_0\" with lines ls 2, "
        gp +=      "'V(r)_push.dat' using (pos($0)):(pot($1)) title \"V_max\" with lines ls 3 \n"

    # t-axis
    if ( do_pot ):
        gp += "reset \n"
        gp += "set style line 1 linetype  1 linewidth 4 linecolor rgb \"blue\" \n"
        gp += "set style line 2 linetype  1 linewidth 4 linecolor rgb \"red\" \n"
        gp += "set size m * blowup, sz \n"
        gp += "set origin (sz-m-d3), 0 \n"
        gp += "set tmargin 6.9 \n"
        gp += "set bmargin 13.90 \n"
        gp += "set lmargin 3.5 \n"
        gp += "set key at screen (sz-0.3*m), screen 0.08*sz \n"
        gp += "unset ytics  \n"
        gp += "unset ylabel \n"
        gp += "set y2tics \n"
        gp += "set y2label 'Zeit in fs' \n"
        gp += "set y2range [0:" + str( Nt * dt * CONV_au_fs ) + "]\n"
        gp += 'set format x "%.1tx10^{%2T}" \n' 
        gp += "set xtics \n"
        gp += "set xlabel 'Emissionsrate in 1/fs' offset 0,-2.5 \n"
        gp += "set x2tics rotate \n"
        gp += "set xtics rotate \n"
        gp += "set x2label 'Potential in eV' offset 0,2 \n"
        tVmin_ = tVmin * dt * CONV_au_fs
        tVmax_ = tVmax * dt * CONV_au_fs
        gp += "set obj 1 circle at second " + str( Vmin * CONV_au_eV ) + "," + str( tVmin_ ) + " size screen 0.008 front fillcolor rgb \"blue\" fillstyle solid 1.0 \n"
        gp += "set obj 2 circle at second " + str( Vmax * CONV_au_eV ) + "," + str( tVmax_ ) + " size screen 0.008 front fillcolor rgb \"dark-magenta\" fillstyle solid 1.0 \n"
        gp += "set obj 3 circle at second " + str( tc_V * CONV_au_eV ) + "," + str( t_charge ) + " size screen 0.008 front fillcolor rgb \"black\" fillstyle solid 1.0 \n"
        gp += "plot 'pulse_shape.dat' using ($1 * " + str( CONV_au_eV ) + "):($0 * " + str( dt * CONV_au_fs ) + ") axis x2y2 title 'Pulse-Potential' with lines ls 1,"
        gp +=       "'current.dat' using ($3 / " + str( CONV_au_fs ) + "):($1 * " + str( CONV_au_fs ) + ") axis x1y2 title 'Tunnel Strom' with lines ls 2\n"

# Palette definition
    #pal_c = [0,0.03125,0.0625,0.09375,0.125,0.15625,0.1875,0.21875,0.25,0.28125,0.3125,0.34375,0.375,0.40625,0.4375,0.46875,0.5,0.53125,0.5625,0.59375,0.625,0.65625,0.6875,0.71875,0.75,0.78125,0.8125,0.84375,0.875,0.90625,0.9375,0.96875,1]
    #pal_r = [0.2314,0.2667,0.302,0.3412,0.3843,0.4235,0.4667,0.5098,0.5529,0.5961,0.6392,0.6824,0.7216,0.7608,0.8,0.8353,0.8667,0.898,0.9255,0.9451,0.9608,0.9686,0.9686,0.9686,0.9569,0.9451,0.9255,0.898,0.8706,0.8353,0.7961,0.7529,0.7059]
    #pal_g = [0.298,0.3529,0.4078,0.4588,0.5098,0.5569,0.6039,0.6471,0.6902,0.7255,0.7608,0.7882,0.8157,0.8353,0.851,0.8588,0.8667,0.8471,0.8275,0.8,0.7686,0.7333,0.6941,0.651,0.6039,0.5529,0.498,0.4392,0.3765,0.3137,0.2431,0.1569,0.0157]
    #pal_b = [0.7529,0.8,0.8431,0.8824,0.9176,0.9451,0.9686,0.9843,0.9961,1,1,0.9922,0.9765,0.9569,0.9333,0.902,0.8667,0.8196,0.7725,0.7255,0.6784,0.6275,0.5804,0.5294,0.4824,0.4353,0.3882,0.3451,0.302,0.2588,0.2196,0.1843,0.149]

    # Matrix
    gp += "reset \n"
    gp += "set size (sz-d1-m), (sz-m-d4) \n"
    gp += "set origin d1, m  \n"
    gp += "unset xtics \n"
    gp += "unset ytics \n"
    gp += "set xrange [0:" + str( Nr ) + "]\n"
    gp += "set yrange [0:" + str( Nt ) + "]\n"
    if ( do_square ):
        gp += "set cbrange [" + str( math.log10( TINY ) ) + ":" + str( math.log10( HUGE ) ) + "] \n"
        gp += "set colorbox horizontal user origin (d1+cbar_hborder),(sz-d4+cbar_vborder) size (sz-m-d1-2*cbar_hborder), (d4-2*cbar_vborder) front \n"
        gp += 'set label "log_{10}( {/Symbol-Oblique Y}^2 nm )" at screen 0.055*sz, screen 0.95*sz front \n'
        #gp += "unset colorbox \n"
        gp += "set pm3d map \n"
        #gp += "set palette defined("
        #for ci in range( len( pal_c ) ):
        #    gp += str( math.log10(TINY) + ci * ( math.log10(HUGE) - math.log10(TINY) ) / ( len(pal_c) - 1 ) )
        #    gp += " " + str( pal_r[ci] ) + " " + str( pal_g[ci] ) + " " + str( pal_b[ci] ) + ", "
        #gp += "\n" 
        #gp += "set palette grey negative\n"
        gp += "set palette defined("
        gp += str( math.log10(TINY) )       + " 1 1 1, "
        gp += str( math.log10(HUGE)-0.03 )  + " 0 0 0, "
        gp += str( math.log10(HUGE) )       + " 1 0.8 0.8) \n"
        gp += "plot 'matrix.dat' matrix with image\n"
    else:
        gp += "plot 'matrix.dat' with image\n"

    proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
    print proc.communicate( gp )[0]

#####################################################################################
#                   Movie frames                                                    #
#####################################################################################

    if ( not do_movie ):
        sys.exit( 0 )

    gp1 = "set term png size 800, 600\n"
    gp1 += "set y2label  \"Amplitude\"\n"
    gp1 += "set yrange [-1.4:0.1]\n"
    gp1 += "set ylabel \"E in Ha\"\n"
    gp1 += "set y2range ["+str( TINY )+":" + str( HUGE ) + "]\n"
    gp1 += "set xlabel \"r in nm\"\n"
    gp1 += "set logscale y2\n"
    gp1 += "#set key off\n"
    gp1 += "#t = sprintf(\"%5.3g\",time)\n"
    gp1 += "set output '"
    
    gp2 = "plot \"-\" using ($0 * " + str(dr) + "):1 axis x1y2 title \"Psi^2\" with lines"
    if ( do_pot ):
        gp2 += ", '' using ($0 * " + str(dr) + "):1 axis x1y1 title \"V_real\" with lines\n"
    else:
        gp2 += "\n"
    
    if ( do_pot ):
        for i in range( Nt ):
            si = str('%06d' % i)
            print si
            stream_in = gp1 + "plot_movie_frame_" + si + ".png' \n" + gp2 
            
            for x in range( Nr ):
                stream_in += str( psi[i][x] ) + "\n"
            stream_in += "e\n"
            for x in range( Nr ):
                stream_in += str( pot[i][x] ) + "\n"
            stream_in += "e\n"
            
            proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
            print proc.communicate( stream_in )[0]
    else:
        for i in range( Nt ):
            si = str('%06d' % i)
            print si
            stream_in = gp1 + "plot_movie_frame_" + si + ".png' \n" + gp2 

            for x in range( 0, Nr ):
                stream_in += str( psi[i][x] ) + "\n"

            proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
            print proc.communicate( stream_in )[0]



    subprocess.call( ["mencoder",    "mf://./plot_movie_frame_*",
                                   "-mf", "type=png:w=800:h=600:fps=10",
                                    "-ovc", "lavc",
                                    "-lavcopts", "vcodec=mpeg4",
                                    "-oac", "copy",
                                    "-o", filename + ".avi" ] )
    subprocess.call( [ "bash", "-c", "rm plot_movie_frame_*" ] )

if __name__ == '__main__':
    main()

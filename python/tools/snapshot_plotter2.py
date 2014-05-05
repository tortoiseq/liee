'''
Created on 04-Feb-2012
forked for simplification version on Sep 25, 2013

Tool for plotting the simulation snapshot files taken by liee::Obs_Snapshot_WF.
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

def multiplot():
    pass

def main():
    # parse arguments, load datafile, find maxima of Psi, minima of V
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:t:r:l", ["file=","tiny=","huge=","logscale","complex"])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    TINY = float('nan')
    HUGE = float('nan')
    filename = ""
    logscale = False
    do_square = True
    t_range = False 
    r_range = False
    
    for o, a in opts:
        if o in ("-f", "--file"):
            filename = a
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
        elif o in ("--complex"):
            do_square = False
        else:
            assert False, "unhandled option"

#####################################################################################
#                       Acquire Data                                                #
#####################################################################################

    r_0 = 0
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
        elif ( line.find( "t_charge\t" ) > -1 ):
            i = line.find( "t_charge\t" )
            t_charge = float( line[( i + len("t_charge\t") ):] ) 
    fobj.close()
    
    if ( do_square == True ):
        psi = numpy.loadtxt( filename, numpy.float, '#', )  # @UndefinedVariable
    else:
        num_col = len( open( filename ).readline().split( "\t" ) )
        psi = numpy.loadtxt( filename, numpy.float, comments='#', usecols = range( 0, num_col-1, 2 ) )  # @UndefinedVariable
        psi_im = numpy.loadtxt( filename, numpy.float, comments='#', usecols = range( 1, num_col-1, 2 ) )  # @UndefinedVariable

    Nt = psi.shape[0]
    Nr = psi[0].shape[0]
    print(Nt)
    print(Nr)
    psi_log = numpy.zeros( (Nt, Nr) )  # @UndefinedVariable
    dr = r_range / ( Nr - 1 )
    dt = t_range / ( Nt - 1 )
    print( "t_charge = " + str(t_charge) )
    print( "t_range = " + str(t_range) )
    print( "Nt = " + str(Nt) )
    tc_i = int( t_charge / CONV_au_fs / dt ) + 1 
    print( "tc_i = " + str(tc_i) )
    
    # auto-scale margins
    mini = float("inf")
    maxi = float("-inf")
    for ti in range( Nt ):
        for ri in range( Nr ):
            absr = abs( psi[ti][ri] )
            mini = absr if ( absr < mini ) else mini
            maxi = absr if ( absr > maxi ) else maxi
            if ( not do_square ):
                psi[ti][ri] = psi[ti][ri] / CONV_au_nm
                absi = abs( psi_im[ti][ri] )
                mini = absi if ( absi < mini ) else mini
                maxi = absi if ( absi > maxi ) else maxi
    TINY = mini if ( math.isnan( TINY ) ) else TINY
    HUGE = maxi if ( math.isnan( HUGE ) ) else HUGE
    print( "Bounds:  [ " + str( TINY ) + " : " + str( HUGE ) + " ]\n")
    
    tVmin = 0
    tVmax = 0
    Vmin = float("inf")
    Vmax = float("-inf")
    lib = ctypes.cdll.LoadLibrary('/usr/local/lib/libliee.so')
    lib.init_potential()
    pot = numpy.zeros( (Nt, Nr) )  # @UndefinedVariable
    c_Vr = ctypes.c_double( 0.0 );
    c_Vi = ctypes.c_double( 0.0 );
    for ti in range( Nt ):
        for ri in range( Nr ):
            c_r = ctypes.c_double( r_0 + ri * dr );
            c_t = ctypes.c_double( ti * dt );
            lib.calc_potential( c_r, c_t, ctypes.byref(c_Vr), ctypes.byref(c_Vi) )
            pot[ti][ri] = c_Vr.value
            if ( pot[ti][ri] > Vmax ):
                Vmax = pot[ti][ri]
                tVmax = ti;
            if ( pot[ti][ri] < Vmin ):
                Vmin = pot[ti][ri]
                tVmin = ti;

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
            numpy.savetxt( "matrix.dat", psi_log )  # @UndefinedVariable
        else:
            numpy.savetxt( "matrix.dat", psi )  # @UndefinedVariable
    else:
        outfile.close()
        
    r_split = 0.0
    i_split = int( ( r_split / CONV_au_nm - r_0 ) / dr )
    print( str(i_split) + " of "+ str(Nr) )
    V_of_t = numpy.zeros( Nt )  # @UndefinedVariable
    V_of_r = numpy.zeros( Nr )  # @UndefinedVariable
    for i in range( Nt ):
        V_of_t[i] = pot[i][i_split]
    numpy.savetxt( "pulse_shape.dat", V_of_t )  # @UndefinedVariable
    
    tc_V = pot[tc_i][i_split]
    for i in range( Nr ):
        V_of_r[i] = pot[tc_i][i]
    numpy.savetxt( "V(t0,r).dat", V_of_r )  # @UndefinedVariable

    for i in range( Nr ):
        V_of_r[i] = pot[tVmin][i]
    numpy.savetxt( "V(r)_pull.dat", V_of_r )  # @UndefinedVariable

    for i in range( Nr ):
        V_of_r[i] = pot[tVmax][i]
    numpy.savetxt( "V(r)_push.dat", V_of_r )  # @UndefinedVariable

#####################################################################################
#                   Potential slices                                                #
#####################################################################################
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
    gp =  "set terminal postscript enhanced eps color solid 'Helvetica' 20 \n"
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
    gp += "set style line 1 linetype  1 linewidth 4 linecolor rgb \"blue\" \n"
    gp += "set style line 2 linetype  1 linewidth 4 linecolor rgb \"black\" \n"
    gp += "set style line 3 linetype  1 linewidth 4 linecolor rgb \"dark-magenta\" \n"
    gp += "set size sz, m * blowup \n"
    gp += "set lmargin 8.2 \n"
    gp += "set rmargin 33 \n"
    gp += "set bmargin 4.0 \n"
    gp += "set origin 0, 0 \n"
    gp += "set ylabel 'Potential in eV' \n"
    gp += "set ytics nomirror\n"
    gp += "set yrange [" + str( (Vmin - 0.03*(Vmax-Vmin)) * CONV_au_eV ) + ":" + str( (Vmax + 0.03*(Vmax-Vmin)) * CONV_au_eV ) + "]\n"
    gp += "set xrange [" + str( r_0 * CONV_au_nm ) + ":" + str( (r_0 + Nr * dr) * CONV_au_nm ) + "]\n"
    gp += "plot 'V(r)_pull.dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_min\" with lines ls 1,"
    gp +=      "'V(r)_push.dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_max\" with lines ls 3,"
    gp +=      "'V(t0,r).dat' using (" + str( r_0 * CONV_au_nm ) + " + $0 * " + str( dr * CONV_au_nm ) + "):($1 * " + str( CONV_au_eV ) + ") title \"V_0\" with lines ls 2\n"

    # t-axis
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
        gp += "set pm3d map \n"
        gp += "set palette defined("
        gp += str( math.log10(TINY) )       + " 1 1 1, "
        gp += str( math.log10(HUGE)-0.03 )  + " 0 0 0, "
        gp += str( math.log10(HUGE) )       + " 1 0.8 0.8) \n"
        gp += "plot 'matrix.dat' matrix with image\n"
    else:
        gp += "plot 'matrix.dat' with image\n"

    proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
    print proc.communicate( gp )[0]

if __name__ == '__main__':
    main()

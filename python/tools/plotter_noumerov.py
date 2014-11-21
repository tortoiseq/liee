'''
Old and busted, not sure if still applicable. Has been unused for a long time.
Created on 04-Feb-2012
Tool for plotting the wave functions from Noumerov debug scans
'''
import numpy
import subprocess
import math
import sys


def main():
    
    TINY = 1e-10
    HUGE = 1e-1
    tend = 0 
    rend = 0
    filename = "wave_function_ini.dat"
#####################################################################################
#                       Acquire Data                                                #
#####################################################################################

    fobj = open("summary.txt", "r")
    for line in fobj: 
        if ( line.find( "r_end\t" ) > -1 ):
            i = line.find( "r_end\t" )
            rend = float( line[( i + len("r_end\t") ):] )
        elif ( line.find( "t_end\t" ) > -1 ):
            i = line.find( "t_end\t" )   
            tend = float( line[( i + len("t_end\t") ):] )
    fobj.close()

    table = numpy.loadtxt( filename, numpy.float, '#', ) 
    Nt = table.shape[0]
    Nr = table[0].shape[0]
    if ( tend == 0 ):
        tend = 1.0
    if ( rend == 0 ):
        rend = 1.0
        
    psi = numpy.zeros( (Nt, Nr) )
    for ti in range( Nt ):
        for ri in range( Nr ):
            psi[ti][ri] = table[ti][ri]
    numpy.savetxt( "matrix.dat", psi )

    #for ti in range( Nt ):
    #    for ri in range( Nr ):
    #        if psi[ti][ri] < TINY:
    #            psi[ti][ri] = math.log10( TINY )
    #        else:
    #            psi[ti][ri] = math.log10( psi[ti][ri] )

    max = -6.66e66
    for i in range( Nt ):
        for x in range( Nr ):
            if ( psi[i][x] > max ):
                max = psi[i][x]
        
    r_scale = rend / ( Nr - 1 )
    #t_scale = tend / ( Nt - 1 )
    
#####################################################################################
#                   Matrix plot                                                     #
#####################################################################################    
    Nr = 1000;
    
    gnuplot0 = "set term png size 1200, 800\n\
#set ylabel \"E (Ha)\"\n\
#set xlabel \"r (a0)\"\n\
#set key off\n\
#set title \"\"\n\
set xrange [0:" + str( Nr ) + "]\n\
set yrange [0:" + str( Nt ) + "]\n\
#set cbrange [" + str( math.log10( TINY ) ) + ":" + str( math.log10( HUGE ) ) + "] \n\
set pm3d map\n\
set output 'matrix.png'\n\
plot \"matrix.dat\" matrix with image\n"

    proc = subprocess.Popen( "gnuplot", stdin=subprocess.PIPE, stdout=subprocess.PIPE )
    print proc.communicate( gnuplot0 )[0]

    #sys.exit( 0 )

#####################################################################################
#                   Movie frames                                                    #
#####################################################################################
    gp1 = "set term png size 800, 600\n"
    gp1 += "set y2label  \"Amplitude\"\n"
    gp1 += "set yrange [-1.4:0.1]\n"
    gp1 += "set ylabel \"E in Ha\"\n"
    gp1 += "#set y2range [0:" + str( max ) + "]\n"
    gp1 += "set xlabel \"r in nm\"\n"
    gp1 += "#set logscale y2\n"
    gp1 += "#set key off\n"
    gp1 += "#t = sprintf(\"%5.3g\",time)\n"
    gp1 += "set output '"
    
    gp2 = "plot \"-\" using ($0 * " + str(r_scale) + "):1 axis x1y2 title \"Psi^2\" with lines \n"
    #gp2 += ", '' using ($0 * " + str(r_scale) + "):1 axis x1y1 title \"V_real\" with lines\n"
    
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

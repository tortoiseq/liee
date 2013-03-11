'''
Created on Nov 28, 2012

@author: quark
'''

import numpy
import sys
import os
import subprocess
import getopt
import math
import ctypes
import re
from GnuPGInterface import GnuPG



def arithmetic_mean( list ):
    return sum( list ) / len( list )

def variance( list ):
    var = 0
    mean = arithmetic_mean( list )
    for x in list:
        var += ( x - mean) ** 2 
    return var

def variance_( mean, list ):
    var = 0
    for x in list:
        var += ( x - mean) ** 2 
    return var

def main():
    pass

    # parse arguments, load datafile, find maxima of Psi, minima of V
    try:
        opts, args = getopt.getopt( sys.argv[1:], "f:d", ["file=","dims="] )  # --dims=1,3,4 [1,3,4]
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    data = []
    filename = "swarming.dat";        
    dims = [];

    for o, a in opts:
        if o in ("-f", "--file"):
            filename = a
        elif o in ("-d", "--dims"):
            dims = re.findall( r'\w+', a )
            for di in range( len( dims ) ):
                dims[di] = int( dims[di] )        
        else:
            assert False, "unhandled option"
    dim = len( dims )  
            
            
    swarm = open( filename, "r" )
    
    # read
    di = -1
    for line in swarm: 
        data.append( [] )
        di += 1
        tupels = line.split("\t")
        ti = -1
        for t in tupels:
            data[di].append( [] )
            ti += 1
            coordinates = t.split("|")
            for dptr in dims:
                if dptr >= len( coordinates ):
                    print( "index error: " + str(dptr) + " > " + str( len(coordinates) ) + "\t@ line " + str(di) + " - " + str(ti) )
                    exit( 1 )
                data[di][ti].append( float( coordinates[dptr] ) )   # only pick the user-selected dimensions
    swarm.close()
    
    avg_range = 3
    zoom_range = 1
    mini = range( len( data ) )
    maxi = range( len( data ) )
    for row_i in range( len( data ) ):
        mini[row_i] = range( dim )
        maxi[row_i] = range( dim )
        for d in range( dim ):
            sample = []
            for row_j in range( max( 0, row_i - avg_range ), min( len( data ), row_i + avg_range + 1 ) ):
                for point in data[row_j]:
                    sample.append( point[d] )
            m = arithmetic_mean( sample )
            v = variance_( m, sample )
            mini[row_i][d] = m - zoom_range * math.sqrt( v )
            maxi[row_i][d] = m + zoom_range * math.sqrt( v )
    
    #interpolation to grid
    gridN = 30
    grid = numpy.zeros( (gridN, gridN) )
    for xi in range( gridN ):
        x = mini[0][0] + xi * ( maxi[0][0] - mini[0][0] ) / ( gridN - 1 ) 
        for yi in range( gridN ):
            y = mini[0][1] + yi * ( maxi[0][1] - mini[0][1] ) / ( gridN - 1 )
            weigth_sum = 0
            z = 0
            bullseye = False
            for row in data:
                for tuple in row:
                    weight = ( x - tuple[0] )**2 + ( y - tuple[1] )**2
                    if weight == 0:
                        z = tuple[2]
                        bullseye = True
                        break
                    else:
                        weight = 1.0 / weight
                        weigth_sum += weight
                        z += weight * tuple[2]
                if bullseye: break
            grid[xi][yi] = z / weigth_sum 
    numpy.savetxt( "grid.dat", grid )
            
    #colormap instead of 3d
    gplot = open( filename+"_.gplot", "w" )
    gnuplot0 = "set terminal gif enhanced size 640,480 animate delay 10\n\
#set ylabel \"E (Ha)\"\n\
#set xlabel \"r (a0)\"\n\
#set key off\n\
#set title \"\"\n\
set pm3d map\n\
set output 'animation_.gif'\n\
plot \"grid.dat\" matrix with image\n"
    gplot.write( gnuplot0 )
    gplot.close()

    
    rot1 = 0.5
    rot2 = 0.3
    intras = 10
    coo = ["x", "y", "z"]
                    
    gplot = open( filename+".gplot", "w" )
    gplot.write( "set terminal gif enhanced size 640,480 animate delay 10\n" )
    gplot.write( "set isosamples 30, 30\n" )
    gplot.write( "set key off\n" )
    gplot.write( "set output 'animation.gif'\n" )

    frame = 0
    for di in range( len( data ) ):
        gplot.write( "set xrange [" + str( mini[di][0] ) + ":" + str( maxi[di][0] ) + "]\n" )
        gplot.write( "set yrange [" + str( mini[di][1] ) + ":" + str( maxi[di][1] ) + "]\n" )
        gplot.write( "set zrange [" + str( mini[di][2] ) + ":" + str( maxi[di][2] ) + "]\n" )
        
        #gplot.write( "reset \nset output frame_" )
        #gplot.write( str(di) )
        #gplot.write( ".png\n")
        #arrorow between all pairs of points
        aid = 1
        for pt1 in range( len( data[di] ) ):
            for pt2 in range( pt1 + 1, len( data[di] ) ):
                    point1 = ""
                    point2 = ""
                    N = min( 3, len( dims ) )
                    for d in range( N ):
                        if d > 0:
                            point1 += ","
                            point2 += ","
                        point1 += str( data[di][pt1][d] )
                        point2 += str( data[di][pt2][d] )
                    gplot.write( "set arrow " + str( aid ) + " from " + point1 + " to " + point2 + " nohead\n" )
                    aid += 1
        gplot.write( "splot 'grid.dat' with dots\n" )
        frame += 1
        
        if ( di < len( data ) - 1 ):
            for i in range( 1, intras ):
                for co in range( 3 ):
                    mi = ( ( intras - i ) * mini[di][co] + i * mini[di+1][co] ) / intras
                    ma = ( ( intras - i ) * maxi[di][co] + i * maxi[di+1][co] ) / intras
                    gplot.write( "set " + coo[co] + "range [" + str( mi ) + ":" + str( ma ) + "]\n" )
                gplot.write( "set view " + str( 45 ) + "," + str( frame * rot2 ) + ";replot\n" )
                frame += 1
                
    #pause at the end
    di = len( data ) - 1
    for i in range( 1, intras ):
        gplot.write( "set view " + str( 45 ) + "," + str( frame * rot2 ) + ";replot\n" )
        frame += 1
             
    gplot.close()

    

if __name__ == '__main__':
    main()
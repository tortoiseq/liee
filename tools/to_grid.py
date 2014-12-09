'''
Interpolates scattered data to project it onto a regular grid, which is useful for plotting.
Limitations: - up to two dimensions of size (M,N)
             - very slow for big grids and lots of input vectors
Output to Stdout as a tab-separated table.
'''

import numpy
import sys
import getopt
import math

def main():
    # parse arguments
    try:
        opts, args = getopt.getopt( sys.argv[1:], "x:y:z:N:M:f:r", ["file="] )
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    #defaults
    x = 1
    y = 2
    z = 3
    N = 100
    M = 100
    ipr = 25     #interpolation range
    fname = "scatter.dat"

    for o, a in opts:
        if o in ("-f", "--file"):
            fname = a
        elif o in ("-x",):
            x = int(a)
        elif o in ("-N",):
            N = int(a)
        elif o in ("-y",):
            y = int(a)
        elif o in ("-M",):
            M = int(a)
        elif o in ("-z",):
            z = int(a)
        elif o in ("-r",):
            ipr = int(a)
        else:
            assert False, "unhandled option"

    xmin = 1e66
    xmax = -1e66
    ymin = xmin
    ymax = xmax

    #get ranges in first pass
    fobj = open( fname, "r" )
    for line in fobj:
        cols = line.split("\t")
        tmp = float( cols[x] )
        if tmp > xmax:
            xmax = tmp
        if tmp < xmin:
            xmin = tmp

        tmp = float( cols[y] )
        if tmp > ymax:
            ymax = tmp
        if tmp < ymin:
            ymin = tmp
    fobj.close()

    dx = (xmax-xmin)/(N-1)
    dy = (ymax-ymin)/(M-1)

    grid = numpy.zeros( (N, M) )
    wsum = numpy.zeros( (N, M) )   # should be integers
    maxweight = 1e9

    #second pass: gather points
    fobj = open( fname, "r" )
    for line in fobj:
        cols = line.split("\t")
        xv = float( cols[x] )
        xi = int( (xv - xmin) / dx )
        yv = float( cols[y] )
        yi = int( (yv - ymin) / dy )
        zv = float( cols[z] )
        # attribute to grid points
        for i in range( 1+xi-ipr, 1+xi+ipr ):
            for j in range( 1+yi-ipr, 1+yi+ipr ):
                if 0 <= i < N and 0 <= j < M:
                    w = 0
                    dist = math.sqrt( ( float(i) - (xv-xmin)/dx) ** 2  +  ( float(j) - (yv-ymin)/dy) ** 2 )
                    if dist < 1.0 / maxweight:
                        w = maxweight
                    else:
                        w = 1.0 / dist
                         
                    grid[i][j] += w * zv
                    wsum[i][j] += w
    fobj.close()

    #third pass: weight and print out
    print("##\txmin=" + str(xmin) +";")
    print("##\tymin=" + str(ymin) +";" )
    print("##\txmax=" + str(xmax) +";" )
    print("##\tymax=" + str(ymax) +";" )
    print("##\tdx=" + str(dx) +";" )
    print("##\tdy=" + str(dy) +";" )
    for i in range( N ):
        line = ""
        for j in range( M ):
            zz = 'NaN'
            if wsum[i][j] > 0:
                zz = grid[i][j] / wsum[i][j]
            line += str(zz) + "\t"
        print( line )



if __name__ == '__main__':
    main()
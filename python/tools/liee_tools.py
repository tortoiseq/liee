#!/usr/bin/python2.7
'''
Provides some command line tools to help with gathering data and plotting the data from experiment-runs. 
'''

import sys
import numpy
#import math
import ctypes
import os

CONV_au_m = 5.2917720859e-11
CONV_au_nm = 5.2917720859e-2
CONV_au_eV = 27.211
CONV_au_s = 2.418884326505e-17
CONV_au_fs = 2.418884326505e-2
CONV_au_V_over_m = 5.1421e11
CONV_au_V = CONV_au_V_over_m * CONV_au_m

precision = 5
lib = ctypes.cdll.LoadLibrary('/usr/local/lib/libliee.so')


'''
Simply transposes matrix data of a given text file and writes the result to an output-file.
(For some reason the used numpy routines work (here) only for python2.7)

Problem     commenting lines disappear
Solution    grep '#' input-file >> output-file.dat
'''
def transpose( infile, outfile ):
    data = numpy.loadtxt( infile, numpy.float, '#', )
    numpy.savetxt( outfile, data.transpose(), fmt="%."+str(precision)+"g" )
    return

def matrix_sums( infile, outfile ):
    data = numpy.loadtxt( infile, numpy.float, '#', )
    fobj = open( outfile, "w" )
    for ti in range( data.shape[0] ):
        sum = 0.0
        for i in range( data.shape[1] ):
            sum += data[ti][i]
        fobj.write( str(sum) + "\n" )
    fobj.close()

def stats( infile, outfile):
    glob_min = 6e66
    glob_max = -6e66
    data = numpy.loadtxt( infile, numpy.float, '#', )
    for y in data.flat:
        if y < glob_min: glob_min = y
        if y > glob_max: glob_max = y
    fobj = open( outfile, "a" )
    fobj.write( "glob_min=" + str( glob_min ) + ";\nglob_max=" + str( glob_max ) + ";\n")
    fobj.close()


def potential( r0, r1, Nr, t0, t1, Nt, outfile, incl_r=0 ):
    if not (incl_r==1 or incl_r==0):
        exit(1)
    r0 = r0 / CONV_au_nm
    r1 = r1 / CONV_au_nm
    t0 = t0 / CONV_au_fs
    t1 = t1 / CONV_au_fs
    if Nr > 1:  dr = (r1 - r0) / (Nr - 1)
    else:       dr = 0
    if Nt > 1:  dt = (t1 - t0) / (Nt - 1)
    else:       dt = 0
    
    lib.init_potential()
    pot = numpy.zeros( (Nr, Nt+incl_r) )
    c_Vr = ctypes.c_double( 0.0 );
    c_Vi = ctypes.c_double( 0.0 );
    V_min = 6e66
    V_min_t = 0
    V_max = -6e66
    V_max_t = 0
    for ti in range( Nt ):
        for ri in range( Nr ):
            c_r = ctypes.c_double( r0 + ri * dr );
            c_t = ctypes.c_double( t0 + ti * dt );
            lib.calc_potential( c_r, c_t, ctypes.byref(c_Vr), ctypes.byref(c_Vi) )
            if incl_r==1 and ti==0:   pot[ri][0] = c_r.value * CONV_au_nm
            pot[ri][ti+incl_r] = c_Vr.value * CONV_au_eV
            if pot[ri][ti+incl_r] < V_min:
                V_min = pot[ri][ti+incl_r]
                V_min_t = ti
            if pot[ri][ti+incl_r] > V_max:
                V_max = pot[ri][ti+incl_r]
                V_max_t = ti
    numpy.savetxt( outfile, pot, fmt="%."+str(precision)+"g" )
    os.system("sync")
    fobj = open( outfile, "a" )
    fobj.write( "##\tt0=" + str( t0 * CONV_au_fs ) + ";\n" )
    fobj.write( "##\tt1=" + str( t1 * CONV_au_fs ) + ";\n" )
    fobj.write( "##\tNt=" + str( Nt ) + ";\n" )
    fobj.write( "##\tr0=" + str( r0 * CONV_au_nm ) + ";\n" )
    fobj.write( "##\tr1=" + str( r1 * CONV_au_nm ) + ";\n" )
    fobj.write( "##\tNr=" + str( Nr ) + ";\n" )
    fobj.write( "##\tVmin=" + str( V_min ) + ";\n" )
    fobj.write( "##\tVmin_i=" + str( V_min_t ) + ";\n" )
    fobj.write( "##\tVmax=" + str( V_max ) + ";\n" )
    fobj.write( "##\tVmax_i=" + str( V_max_t ) + ";\n" )
    fobj.close()

def plot( module_id ):
    if ( module_id == "all" ):
        module_id = -1
    plot_id = ctypes.c_int( int(module_id) )
    lib.plot( plot_id )


if __name__ == '__main__':
    if sys.argv[1] == "transpose":
        transpose( sys.argv[2], sys.argv[3] )

    elif sys.argv[1] == "stats":
        if ( len(sys.argv) > 4 and sys.argv[4] == "skipcol0"):
            os.system("sed -e \"s/[\t[:space:]]\+/\t/g\" '"+ sys.argv[2] +"' | cut -f 2- - >'" + sys.argv[2] + ".skipcol0'")
            stats( sys.argv[2]+".skipcol0", sys.argv[3] )
        else:
            stats( sys.argv[2], sys.argv[3] )

    elif sys.argv[1] == "potential":
        if len(sys.argv) > 9:
            if sys.argv[9].upper()=="TRUE":
                potential( float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7]), sys.argv[8], 1 )
            elif sys.argv[9].upper()=='FALSE':
                potential( float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7]), sys.argv[8], 0 )
        else:
            potential( float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7]), sys.argv[8], 0 )

    elif sys.argv[1] == "plot":
        plot( sys.argv[2] )

    elif sys.argv[1] == "sums":
        matrix_sums( sys.argv[2], sys.argv[3] )

    else:
        print("usage:\n" +
              "* liee_tools transpose \t input-file output-file \n" +
              "* liee_tools stats \t input-file output-file \n" +
              "* liee_tools stats \t input-file output-file skipcol0\n" +
              "* liee_tools potential \t r0 r1 Nr t0 t1 Nt outfile [true/false]\n" +
              "* liee_tools plot \t [all|module_serial]\n" +
              "* liee_tools sums \t input-file output-file \n" +
              "" )

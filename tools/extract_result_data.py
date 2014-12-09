'''
Extracts specific values from a set of finished results.

An LIEE-experiment typically runs many instances with modified parameters resulting in a bunch of directories
in the servers upload-path containing the individual data-files. The input parameters and some results are listed
in the "summary.txt" for each run. "extract_result_data" currently can only access those values listed in the summary file.

Optionally predefined plotting routines can be executed within each result-directory. 

The output is a table with a row for each experiment-run and the unnamed parameters in the given order ("params"-option)
'''

import numpy
import sys
import os
import subprocess
import getopt
import math
import ctypes
import re

def main():
    # parse arguments, load datafile, find maxima of Psi, minima of V
    try:
        opts, args = getopt.getopt( sys.argv[1:], "e:p", ["experiment=","params=","plot-each"] )
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    plot_each = False
    ex = ""  # experiment name
    ps = []  # list of parameters to extract
    dirs = []  # remember appropriate directories for the plotting
    fi = -1  # file index
    data = []

    for o, a in opts:
        if o in ("-e", "--experiment"):
            ex = a
        elif o in ("-p", "--params"):
            ps = re.findall( r'\w+', a )
            print( ps )
        elif o in ("--plot-each",):
            plot_each = True        
        else:
            assert False, "unhandled option"


    for dirname, dirnames, filenames in os.walk( '.' ):
    # filter directories according to experiment
        if ( dirname.find( ex ) == -1 ):
            continue

        for filename in filenames:
            if filename == "summary.txt":
                fi += 1
                fobj = open( os.path.join( dirname, filename ), "r" )
                data.append( range( len( ps ) ) )
                for di in range( len( ps ) ):
                    data[fi][di] = float('nan')
                    
                for line in fobj:
                    for pi in range( len( ps ) ): 
                        key = " " + ps[pi] + "\t"
                        if ( line.find( key ) > -1 ):
                            i = line.find( key )
                            data[fi][pi] = float( line[( i + len( key ) ):] )
                fobj.close()
                
                if plot_each:
                    dirs.append( dirname )

    # create the plots specified in liee_parameter.xml 
    if plot_each:
        for dir in dirs:
            pass # nada
            sufix = "_"
            i = dir.find( "_" )
            if i > -1:
                sufix = dir[ (i+1): ]
            print( dir )
            os.chdir( dir )
            #TODO use new plotting tools
            os.chdir( ".." )


    # OUTPUT: print to stdout
    for fi in range( len( data ) ):
        line = str( fi )
        for x in data[fi]:
            line += "\t" + str( x )
        print( line )



if __name__ == '__main__':
    main()
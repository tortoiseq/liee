'''
Created on Nov 26, 2012

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

def main():
    # parse arguments, load datafile, find maxima of Psi, minima of V
    try:
        opts, args = getopt.getopt( sys.argv[1:], "e:p", ["experiment=","params=","plot-each","tiny=","huge="] )
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    plot_each = False
    ex = ""
    ps = []
    dirs = []
    fi = -1
    data = []
    tiny = "1e-9"
    huge = "1e-1"
        
    for o, a in opts:
        if o in ("-e", "--experiment"):
            ex = a
        elif o in ("-p", "--params"):
            ps = re.findall( r'\w+', a )
            print( ps )
        elif o in ("--plot-each",):
            plot_each = True        
        elif o in ("--tiny",):
            tiny = a        
        elif o in ("--huge",):
            huge = a        
        else:
            assert False, "unhandled option"
    
    
    for dirname, dirnames, filenames in os.walk( '.' ):
    #filter directories according to experiment
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
                
    # run snapshot plotter
    if plot_each:
        for dir in dirs:
            sufix = "_"
            i = dir.find( "_" )
            if i > -1:
                sufix = dir[ (i+1): ]
            print( dir )
            os.chdir( dir )
            subprocess.call( [ "python", "/var/www/boinc/liee/bin/snapshot_plotter.py", "-f", "wf_evolution.dat",
                               "--logscale", "--potential", 
                               "--tiny="+tiny, "--huge="+huge ] )
            #subprocess.call( [ "convert", "multiplot.ps", "multiplot" + sufix + ".png" ] )
            os.chdir( ".." )
                

    #print stdout
    for fi in range( len( data ) ):
        line = str( fi )
        for x in data[fi]:
            line += "\t" + str( x )
        print( line )



if __name__ == '__main__':
    main()
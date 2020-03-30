#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__           = "Dilawar Singh"
__copyright__        = "Copyright 2017-, Dilawar Singh"
__version__          = "1.0.0"
__maintainer__       = "Dilawar Singh"
__email__            = "dilawars@ncbs.res.in"
__status__           = "Development"

import sys
import os
import pyabf
import numpy as np
import pandas as pd
try:
    import pathlib
except ImportError as e:
    print( 'python > 3.4 is required' )
    quit(1)

args_ = None

def extractData(abf):
    df = pd.DataFrame()
    header = abf.headerText
    for sl in abf.sweepList:
        abf.setSweep(sl)
        x, y = abf.sweepX, abf.sweepY
        df[ 'Time' ] = x
        df[ 'Trace %s'%sl ] = y
    return df

def plotFigure( df, outfile ):
    import matplotlib.pyplot as plt
    df.plot()
    plt.savefig( outfile )
    print( '[INFO] Saved to %s' % outfile )

def saveData( df, outfile ):
    odir = pathlib.Path(outfile).parent
    odir.mkdir(parents=True, exist_ok=True)
    df.to_csv( outfile, index = False )
    print( '[INFO] Saved to %s' % outfile )

def find_files( directory ):
    res = []
    for d, sd, fs in os.walk(directory):
        for f in fs:
            ext = f.split( '.' )[-1].lower()
            if ext == 'abf':
                res.append(os.path.join(d,f))
    return res

def relative_path( filepath, relativeTo = None ):
    if not relativeTo:
        return filepath
    return os.path.realpath(filepath).lstrip(os.path.realpath(relativeTo))

def process_files( files ):
    global args_
    if args_.dir:
        resdir = os.path.join(args_.dir, '_CSV')
        if not os.path.isdir(resdir):
            os.makedirs(resdir)
    else:
        resdir = os.path.dirname(args_.input)

    for f in files:
        abf = pyabf.ABF( f )
        df = extractData(abf)
        saveData(df, os.path.join(resdir
            , '%s.csv' % relative_path(f, args_.dir)))

def main(args):
    global args_
    args_ = args
    if args_.dir:
        files = find_files( args_.dir )
        print( '[INFO] Found %d files' % len(files))
        process_files(files)
    else:
        process_files( [ args_.input ] )


if __name__ == '__main__':
    import argparse
    # Argument parser.
    description = '''ABF format to other formats.'''
    parser = argparse.ArgumentParser(description=description)
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--input', '-i'
        , required = False, help = 'Input file'
        )
    g.add_argument( '--dir', '-d'
            , required = False, help = 'Directory of data'
            )
    parser.add_argument('--output', '-o'
        , required = False, help = 'Output file. Valid only when -i is used.'
        )
    class Args: pass
    args = Args()
    parser.parse_args(namespace=args)
    main(args)
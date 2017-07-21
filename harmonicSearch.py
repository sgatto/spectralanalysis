#!/usr/bin/python

###loading shell commands
import os, os.path, sys, shutil
from scipy.fftpack import fftn
# from PyQt4.QtGui import *
import struct
import numpy as np
import math
from test import *
from spaceAnalysisOnly import *
#import sys
import textwrap
import argparse
# ... do something with args.output ...
    # ... do something with args.verbose ..
#from evtk.hl import gridToVTK

def run(args):
    myxmin=args.myxmin
    myxmax=args.myxmax
    myymin=args.myymin
    myymax=args.myymax
    mycomp=args.mycomp
    print str(args)
    #print 'Number of arguments:', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)

    myname = args.myname
    myanalysis = FieldAnalysis(filename=myname)
    myanalysis.do_fft(zposition=0, comp=mycomp)
    myanalysis.saveffttxt(varname=myname, kxmin=-40, kxmax=0, kymin=-10, kymax=40)

    myanalysis.createCone(nmax=40,phimin=5.5, phimax=10.5, nbin=400)
    myanalysis.analiseConeAndPrint(varname=myname, kxmin=-40, kxmax=0, kymin=-10, kymax=40)
    myanalysis.printConeAnalysis(varname=myname)

    myanalysis.createCone(nmax=40,phimin=2, phimax=12, nbin=400)
    myanalysis.analiseConeAndPrint(varname=myname, kxmin=-40, kxmax=0, kymin=-10, kymax=40)
    myanalysis.printConeAnalysis(varname=myname)




parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        FFT of my files
        --------------------------------
        loop on components (0:x 1:y 2:z)

        '''))
parser.add_argument('-i', type=str, dest='myname', required=True,
                    help='filename of the file to be analised')
parser.add_argument('-xmin', type=float, dest='myxmin', default=-1e5,
                    help='xmin of the box to be analysed')
parser.add_argument('-xmax', type=float, dest='myxmax', default=1e5,
                    help='xmax of the box to be analysed')
parser.add_argument('-ymin', type=float, dest='myymin', default=-1e5,
                    help='ymin of the box to be analysed')
parser.add_argument('-ymax', type=float, dest='myymax', default=1e5,
                    help='ymax of the box to be analysed')
parser.add_argument('-comp', type=int, dest='mycomp', required=True, help='component to analize')
args = parser.parse_args()
run(args)

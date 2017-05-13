#!/usr/bin/python

###loading shell commands
import os, os.path, sys, shutil
from scipy.fftpack import fftn
# from PyQt4.QtGui import *
import struct
import numpy as np
import math
from test import *
from timeanalysisreal import *
#import sys
import textwrap
import argparse
# ... do something with args.output ...
    # ... do something with args.verbose ..
#from evtk.hl import gridToVTK

def run(args):
    mystart=args.mytstart
    mystop=args.mytstop+1
    mycompstart=args.mycstart
    mycompstop=args.mycstop+1
    if(mystop<0):
        mystop=mystart+5
    if(mycompstop<0):
        mycompstop=mycompstart+1
    print str(args)
    #print 'Number of arguments:', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)

    deltaT = 20

    mynames    = ("B_FIELD_", "E_FIELD_", "E_FIELD_","DENS_eleBulk_", "DENS_proBulk_")
    subs       = (	    5   ,      5    ,    	5   ,      5        ,      1         )
    components = (	    2   ,      0    ,    	1   ,      0        ,      0         )
    varnames   = (    "Bz"  ,    "Ex"   ,     "Ey"  ,     "ne"      ,     "np"       )
    #myname = "DENS_eleBulk_"
    #myname = "B_FIELD_"
    #component = 2

    for contatore in range (mycompstart,mycompstop,1):
        myname = mynames[contatore]
        sub=subs[contatore]
        component = components[contatore]
        varname   = varnames[contatore]

        name = ("%s-energy-evolution.txt" % (varname))
        myfile = open(name, 'a')
        if (mystart<=0):
            myfile.write("#1time,   2Etot,  3EOmega2\n")

        for start in range(mystart, mystop, 50):
            basetime = start-deltaT*0.5
            if(basetime<0):
                basetime = 0
            endtime = basetime + deltaT
            centertime = 0.5*(basetime + endtime)
            myanalysis = FieldAnalysis(basetime, endtime, substeps=sub, base=myname, end=".bin.000")
            myanalysis.collect_data()
            myanalysis.setWindow()
            myanalysis.do_fft(zposition=0, comp=component)
            myanalysis.saveVTK3Dfft(varname)
            myanalysis.saveKxOmega(varname)
            if(sub>1):
                myanalysis.saveKxKyatomega(2,varname)
                myanalysis.saveKxKyatomega(1,varname)
            myanalysis.saveKxKyatomega(0,varname)
            myanalysis.printEnergy(zposition=0, comp=component)
            energyAtOmega1=0
            energyAtOmega2=0
            if(sub>1):
                energyAtOmega1=getEnergyAtOmega(1)
                energyAtOmega2=getEnergyAtOmega(2)
                if(component==2):
                    myanalysis.setAllToZeroExceptOmega(2)
                else:
                    myanalysis.setAllToZeroExceptOmega(1)
            myanalysis.printEnergy(zposition=0, comp=component)
            myfile.write("%e, %e, %e, %e\n" % (centertime, myanalysis.totalEnergyFunction, myanalysis.totalEnergyFFT, energyAtOmega1,energyAtOmega2)
        myfile.close()


        #myanalysis.saveKxKyatomega(2.2)
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        FFT of my files
        --------------------------------
        loop on components (0:Bz 1:Ex 2:Ey 3:ne 4:np)
        loop on time windows (delta t=20)
        time increments of 50

        '''))
parser.add_argument('-tstart', type=int, dest='mytstart', required=True,
                    help='the first time of the analysis')
parser.add_argument('-tstop', type=int, dest='mytstop', default=-10,
                    help='the last time to analize')
parser.add_argument('-cstart', type=int, dest='mycstart', required=True, help='integer first component to analize')
parser.add_argument('-cstop', type=int, dest='mycstop', default=-10, help='integer last component to analize')
args = parser.parse_args()
run(args)

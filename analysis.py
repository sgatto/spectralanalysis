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
#from evtk.hl import gridToVTK

def run():
    #print 'Number of arguments:', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    
    deltaT = 20
    
    mynames    = ("B_FIELD_", "E_FIELD_",	"DENS_eleBulk_", "DENS_proBulk_")
    subs       = (	     5   ,    	5    ,	      5        ,      1         )
    components = (	     2   ,    	0    ,	      0        ,      0         )
    varnames   = (    "fftBz",    "fftEx",	        "fftne",         "fftnp")
    #myname = "DENS_eleBulk_"                                                                             
    #myname = "B_FIELD_"                                                                                  
    #component = 2                                                                                        

    for contatore in range (0,1,1):
        myname = mynames[contatore]
        sub=subs[contatore]
        component = components[contatore]
        varname   = varnames[contatore]

        name = ("energy-evolution-%s.txt" % (myname))
        myfile = open(name, 'w')
        myfile.write("#1time,   2Etot,  3EOmega2\n")
        
        for start in range(140, 280, 50):
            basetime = start
            endtime = basetime + deltaT
            centertime = 0.5*(basetime + endtime)
            myanalysis = FieldAnalysis(basetime, endtime, substeps=sub, base=myname, end=".bin.000")
            myanalysis.collect_data()
            myanalysis.setWindow()
            myanalysis.do_fft(zposition=0, comp=component)
            myanalysis.saveVTK3Dfft(varname)
            myanalysis.saveKxOmega()
            if(sub>1):
                myanalysis.saveKxKyatomega(2)
                myanalysis.saveKxKyatomega(1)
            myanalysis.saveKxKyatomega(0)
            myanalysis.printEnergy(zposition=0, comp=component)
            if(sub>1):
                myanalysis.setAllToZeroExceptOmega(2)
            #myanalysis.do_inversefft()
            #myanalysis.saveNewData(0.5*deltaT)
            #myanalysis.saveVTK3Dfft()
            #myanalysis.saveVTK3Dfft(varname, appendix="-mio")
            myanalysis.printEnergy(zposition=0, comp=component)
            myfile.write("%e, %e, %e\n" % (centertime, myanalysis.totalEnergyFunction, myanalysis.totalEnergyFFT))
        myfile.close()


        #myanalysis.saveKxKyatomega(2.2)
run()


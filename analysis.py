#!/usr/bin/python

###loading shell commands
import os, os.path, sys, shutil
from scipy.fftpack import fftn
# from PyQt4.QtGui import *
import struct
import numpy as np
import math
from test import *
from timeanalysis import *
#from evtk.hl import gridToVTK

def run():
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)
    #path = os.getcwd()
    # f        = open(os.path.join(path,'E_FIELD_000.000.bin.000'),'rb')
    component = 0
    basetime = 0
    deltaT = 20
    sub = 5
    endtime = basetime + deltaT
    myname = "E_FIELD_"
    
    for start in range(0, 350, 10):
        print "analizzo il numero %d\n" % start
        basetime = start
        endtime = basetime + deltaT
        myanalysis = FieldAnalysis(basetime, endtime, substeps=sub, base=myname, end=".bin.000", reflect=False)
        myanalysis.collect_data()
        myanalysis.setWindow()
        myanalysis.do_fft(zposition=0, comp=component)
        myanalysis.saveVTK3Dfft()
        myanalysis.saveKxOmega()
run()


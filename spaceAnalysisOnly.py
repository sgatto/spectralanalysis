###!/usr/bin/python
##loading shell commands
import os
import os.path
import glob
import sys
import shutil
import time
from scipy import fftpack
#import scipy.fftpack
# from PyQt4.QtGui import *
import struct
# from scipy import *
import numpy as np
import math
# from evtk.hl import gridToVTK

# from matplotlib import *
# from pylab import *
# import matplotlib as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# ---


def time_filter(loctime, tmax):
    """add envelope time"""
    return 0.5*(1-math.cos(2*np.pi*loctime/tmax))

class mygrid:
    def __init__(self):
        self.init=True
        self.Nx=None
        self.Ny=None
        self.Nz=None
        self.Nt=None
        self.totNr=None
        self.totNk=None

        self.Lx=1
        self.Ly=1
        self.Lz=1
        self.Lt=1
        self.x=None
        self.y=None
        self.z=None
        self.t=None
        self.dx=1
        self.dy=1
        self.dz=1
        self.dt=1

        self.Nkt=None
        self.Nkx=None
        self.Nky=None
        self.Lkx=1
        self.Lky=1
        self.Lkz=1
        self.Lkt=1
        self.kx=None
        self.ky=None
        self.kz=None
        self.kt=None
        self.dkx=1
        self.dky=1
        self.dkz=1
        self.dkt=1

class FieldAnalysis:
    """ la mia classe"""

    def __init__(self, filename):
        self.grid=mygrid()
        self.filename = filename

        self.totalEnergyFunction = 0
        self.totalEnergyFFT = 0

        self.init_values()
        self.print_parameters()
        self.set_frequency()
        self.print_frequency()
        self.collect_data()

    def collect_data(self):
        self.alldata = np.zeros((self.grid.Nz, self.grid.Ny, self.grid.Nx, self.Nc))
        self.alldata = self.analize_field(self.filename)

    def analize_field(self, filename, analize=False):
        f = open(filename, 'rb')
        # - endianness 0:small - 1:big -#

        endianness = struct.unpack('i', f.read(4))[0]

        # - global grid dim -#
        nx = struct.unpack('i', f.read(4))[0]
        ny = struct.unpack('i', f.read(4))[0]
        nz = struct.unpack('i', f.read(4))[0]
        self.grid.Nx = nx
        self.grid.Ny = ny
        self.grid.Nz = nz
        ntot = nx * ny * nz
        # - processor grid -#
        npx = struct.unpack('i', f.read(4))[0]
        npy = struct.unpack('i', f.read(4))[0]
        npz = struct.unpack('i', f.read(4))[0]

        nproc = npx * npy * npz

        # - field components -#
        nc = struct.unpack('i', f.read(4))[0]
        self.grid.Nx = nx
        self.grid.Ny = ny
        self.grid.Nz = nz
        self.Nc = nc

        print nx, ny, nz, nc

        # - grid -> X -#
        # x = np.zeros((Nx))
        # for i in range(0,Nx):
        x = struct.unpack('f' * nx, f.read(4 * nx))

        # - grid -> Y -#
        y = struct.unpack('f' * ny, f.read(4 * ny))
        self.x = x
        # - grid -> Z -#
        z = struct.unpack('f' * nz, f.read(4 * nz))
        self.x = x
        self.y = y
        self.z = z
        # - loop on processors -#
        F = np.zeros((nz, ny, nx, nc))

        print "nprocs=", nproc
        counter = 0
        prog = 0.0
        if(analize):
            return 0
        for nprocessor in range(0, nproc):
            print "proc = ", nprocessor
            # -processor dims -#
            i0 = struct.unpack('i', f.read(4))[0]
            j0 = struct.unpack('i', f.read(4))[0]
            k0 = struct.unpack('i', f.read(4))[0]
            li0 = struct.unpack('i', f.read(4))[0]
            lj0 = struct.unpack('i', f.read(4))[0]
            lk0 = struct.unpack('i', f.read(4))[0]
            # print '>>> ',i0,j0,k0,li0,lj0,lk0
            NN = li0 * lj0 * lk0 * nc
            array = np.array(struct.unpack('f' * NN, f.read(4 * NN))).reshape(lk0, lj0, li0, nc)

            for k in range(0, lk0):
                for j in range(0, lj0):
                    for i in range(0, li0):
                        for c in range(0, nc):
                            F[k + k0, j + j0, i + i0, c] = array[k, j, i, c]
                            counter += li0
                            prog = counter * (100.0 / ntot)

        # np.savetxt( nameOutFile ,F[0,:,:,component],fmt='%15.14e')
        f.close()
        # print "done"
        return F

    def init_values(self):
        self.analize_field(self.filename, analize=True)

        if self.grid.Nx > 1:
            self.grid.dx = self.x[1] - self.x[0]
        else:
            self.grid.dx = 0
        if self.grid.Ny > 1:
            self.grid.dy = self.y[1] - self.y[0]
        else:
            self.grid.dy = 0
        if self.grid.Nz > 1:
            self.grid.dz = self.z[1] - self.z[0]
        else:
            self.grid.dz = 0

        self.grid.Lx = self.grid.dx * self.grid.Nx
        self.grid.Ly = self.grid.dy * self.grid.Ny
        self.grid.Lz = self.grid.dz * self.grid.Nz
        self.grid.totNr = self.grid.Nx * self.grid.Ny * self.grid.Nz

    def print_parameters(self):
        print ("Start Analysis:")
        print ("filename = %s" %self.filename)
        print ("SIZE: [ Lx, Ly, Lz ] = [ %.3f, %.3f, %.3f ]" % (self.grid.Lx, self.grid.Ly, self.grid.Lz))
        print (" Np : [ Nx, Ny, Nz ] = [ %d, %d, %d ]" % (self.grid.Nx, self.grid.Ny, self.grid.Nz))

    def set_frequency(self):

        self.grid.kx = np.fft.fftfreq(self.grid.Nx,d=self.grid.dx)
        self.grid.ky = np.fft.fftfreq(self.grid.Ny,d=self.grid.dy)

        self.grid.kx = np.fft.fftshift(self.grid.kx)
        self.grid.ky = np.fft.fftshift(self.grid.ky)
        #self.grid.kt = np.fft.fftshift(self.grid.kt)

        self.grid.Nkx = self.grid.kx.size
        self.grid.Nky = self.grid.ky.size
        self.grid.totNk = self.grid.Nkx * self.grid.Nky

        self.grid.dkx = 1 / self.grid.Lx
        self.grid.dky = 1 / self.grid.Ly

        self.grid.Lkx = self.grid.dkx * (self.grid.Nx / 2)
        self.grid.Lky = self.grid.dky * (self.grid.Ny / 2)

        self.grid.kz = 0
        self.grid.dkz = 0
        self.grid.Lkz = self.grid.dkz * (self.grid.Nz / 2)


    def print_frequency(self):
        print ("SIZE: [ Lkx, Lky, Lkz ] = [ %6.3f, %6.3f, %6.3f ]" % (self.grid.Lkx, self.grid.Lky, self.grid.Lkz))
        print (" dk : [ dkx, dky, dkz ] = [ %6.3f, %6.3f, %6.3f ]" % (self.grid.dkx, self.grid.dky, self.grid.dkz))

    def do_fft(self, zposition, comp):
        print ("ready for fft3D...")
        dataselect = self.alldata[zposition, :, :, comp]
        self.trasf3D = np.fft.fftn(dataselect, axes=(0,1))
#        self.trasf3D = self.trasf3D/math.sqrt(self.grid.Nx*self.grid.Ny)

        self.shiftedTrasf3D = np.fft.fftshift(self.trasf3D, axes=(0,1))

        print ("DONE fft3D")

    def do_inversefft(self):
        print ("ready for ifft3D...")

        self.mynewData = np.fft.ifftn(self.trasf3D, axes=(1,2))
        self.mynewData = self.mynewData*math.sqrt(self.grid.Nx*self.grid.Ny)
        print ("DONE ifft3D")

    def saveffttxt(self, varname, kxmin=-40, kxmax=0, kymin=-10, kymax=40):
        name = ("%s-2D-fft.txt" % (varname))

        print ("ready for printin on %s..." %(name))
        ikmin = np.maximum((int)((kxmin - self.grid.kx[0])/self.grid.dkx), 0)
        ikmax = np.minimum((int)((kxmax - self.grid.kx[0])/self.grid.dkx), self.grid.Nkx)
        jkmin = np.maximum((int)((kymin - self.grid.ky[0])/self.grid.dky), 0)
        jkmax = np.minimum((int)((kymax - self.grid.ky[0])/self.grid.dky), self.grid.Nky)

        print "self.grid.Nkx = ", self.grid.Nkx, "   self.grid.Nky = ", self.grid.Nky
        print "ikmin = ", ikmin, "   ikmax = ", ikmax, "   jkmin = ", jkmin, "   jkmax = ", jkmax
        f1 = open(name, 'w')
        for j in range(jkmin, jkmax):
            for i in range(ikmin, ikmax):
                f1.write("%.3e, %.3e, %.3e\n" % (self.grid.kx[i], self.grid.ky[j], (np.absolute(self.shiftedTrasf3D[j, i]))**2 ))
        # np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
        f1.close()
        print ("DONE for printin on %s..." %(name))

    def saveNewData(self, varname):

        name = ("%s-x-y.txt" % (varname))
        f1 = open(name, 'w')
        for j in range(0, self.grid.Ny):
            for i in range(0, self.grid.Nx):
                f1.write("%.3e, %.3e, %.3e\n" % (self.x[i], self.y[j], (self.mynewData[itime,j, i]) ) )
            f1.write("\n")
        f1.close()

    def getEnergyAtOmega(self, omegain):

        ifreqin = int(round(omegain/self.grid.dkt))
        print ">>>>>>>ifreqin=" + str(ifreqin) + " omegain=" + str(omegain)
        omegaout = self.grid.kt[ifreqin]
        print " omegaout=" + str(omegaout)
        factor = 8.0/(3.0*self.grid.Nt)*(self.grid.dx * self.grid.dy*2) / (8*math.pi)
        energyAtOmega=0
        for t in range(0, self.grid.Nkt/2):
            if abs(t-ifreqin)<=1:
                energyAtOmega += np.tensordot(self.trasf3D[t,:, :],self.trasf3D[t,:, :].conjugate(),axes=2)
                if t > 0:
                    energyAtOmega += np.tensordot(self.trasf3D[self.grid.Nkt-t,:, :],self.trasf3D[self.grid.Nkt-t,:, :].conjugate(),axes=2)

        return np.absolute(energyAtOmega)*factor

    def createCone(self,nmax,phimin, phimax, nbin):

        print ("creating conditions for the cone...")

        self.nbin = nbin
        self.plotCone = np.zeros(self.nbin)
        self.nmax = nmax
        self.nmin = 0
        self.dkn = nmax*1.0/self.nbin
        self.kns = np.arange(0, self.nmax, self.dkn)
        self.gradphimin = phimin
        self.gradphimax = phimax
        self.phimin = (phimin + 90.0)/180.*np.pi
        self.phimax = (phimax + 90.0)/180.*np.pi
        print "phimin = ", phimin, "   phimax = ", phimax, "   nmax = ", nmax
        print "self.phimin = ", self.phimin, "   self.phimax = ", self.phimax, "   self.nmax = ", self.nmax


    def analiseCone(self):

        print ("start analysis of the cone...")
        print "self.phimin = ", self.phimin, "   self.phimax = ", self.phimax, "   self.nmax = ", self.nmax

        self.plotCone[:] = 0
        for j in range(0, self.grid.Nky):
            ky = self.grid.ky[j]
            for i in range(0, self.grid.Nkx):
                kx = self.grid.kx[i]

                phi = np.arctan2(ky, kx)
                if self.phimin <= phi <= self.phimax:
                    kr = math.sqrt(kx*kx + ky*ky)
                    ikn = int(kr/self.dkn + 0.5)
                    if ikn < self.nbin:
                        self.plotCone[ikn] += ((np.absolute(self.shiftedTrasf3D[j, i]))**2)
        print ("DONE")


    def analiseConeAndPrint(self, varname, kxmin=-40, kxmax=0, kymin=-10, kymax=40):

        print ("start analysis of the cone and print new 2D file...")

        ikmin = np.maximum((int)((kxmin - self.grid.kx[0])/self.grid.dkx), 0)
        ikmax = np.minimum((int)((kxmax - self.grid.kx[0])/self.grid.dkx), self.grid.Nkx)
        jkmin = np.maximum((int)((kymin - self.grid.ky[0])/self.grid.dky), 0)
        jkmax = np.minimum((int)((kymax - self.grid.ky[0])/self.grid.dky), self.grid.Nky)

        print "self.grid.Nkx = ", self.grid.Nkx, "   self.grid.Nky = ", self.grid.Nky
        print "ikmin = ", ikmin, "   ikmax = ", ikmax, "   jkmin = ", jkmin, "   jkmax = ", jkmax
        print "self.phimin = ", self.phimin, "   self.phimax = ", self.phimax, "   self.nmax = ", self.nmax

        name = ("%s-%.1f-%.1f-2DFFT.txt" % (varname, self.gradphimin, self.gradphimax))
        f1 = open(name, 'w')

        self.plotCone[:] = 0

        for j in range(jkmin, jkmax):
            ky = self.grid.ky[j]
            for i in range(ikmin, ikmax):
                kx = self.grid.kx[i]

                phi = np.arctan2(ky, kx)
                if (phi >= self.phimin) and (phi <= self.phimax):
                    kr = math.sqrt(kx*kx + ky*ky)
                    ikn = int(kr/self.dkn + 0.5)
                    if ikn < self.nbin:
                        self.plotCone[ikn] += ((np.absolute(self.shiftedTrasf3D[j, i]))**2)
                    f1.write("%.3e, %.3e, %.3e\n" % (kx, ky, (np.absolute(self.shiftedTrasf3D[j, i]))**2 ) )
                else:
                    f1.write("%.3e, %.3e, %.3e\n" % (kx, ky, 0 ))
        print ("DONE")


    def printConeAnalysis(self, varname):

        print ("print analysis...")
        name = ("%s-%.1f-%.1f.txt" % (varname, self.gradphimin, self.gradphimax))
        f1 = open(name, 'w')

        for i in range(0,self.nbin):
            f1.write("%.3e, %.3e\n" % (self.kns[i], self.plotCone[i]) )
        print ("done")

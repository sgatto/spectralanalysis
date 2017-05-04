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

    def __init__(self, tmin, tmax, substeps, base, end):
        self.grid=mygrid()
        self.base = base
        self.end = end
        self.tmin = int(tmin)
        self.tmax = int(tmax)
        self.substeps = int(substeps)
        if 1000 % substeps:
            print "ERROR: wrong number of timesteps!"
            sys.exit()
        self.grid.Nt = substeps * (self.tmax - self.tmin)
        self.filenames = list()
        self.grid.t = np.zeros((self.grid.Nt))
        self.basename = ("%s%03d" % (self.base, ((self.tmin+self.tmax)/2) ))

        self.totalEnergyFunction = 0
        self.totalEnergyFFT = 0

        self.createFileList()

        self.init_values()
        self.print_parameters()
        self.set_frequency()
        self.print_frequency()

    def analize_field(self, filename):
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
        counter = 0
        prog = 0.0
        for nprocessor in range(0, nproc):
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


    def createFileList(self):
        base = self.base
        end = self.end
        tindex = 0
        for i in range(self.tmin, self.tmax):
            for v in range(0, self.substeps):
                t = i + v * 1.0 / self.substeps

                self.grid.t[tindex] = t
                name = base + ("%03d.%03d" % (i, v * 1000 / self.substeps)) + end
                self.filenames.append(name)
                tindex += 1

        if (self.grid.Nt != len(self.filenames)):
            print ("error self.grid.Nt=%d len(self.filenames))=%d" %(self.grid.Nt, len(self.filenames)))

    def collect_data(self):
        self.alldata = np.zeros((self.grid.Nt, self.grid.Nz, self.grid.Ny, self.grid.Nx, self.Nc))
        for i in range(0, self.grid.Nt):
            # print self.filenames[i]
            F = self.analize_field(self.filenames[i])
            self.alldata[i, :, :, :, :] = F[:, :, :, :]

    def setWindow(self):
        for i in range(0, self.grid.Nt):
            self.alldata[i, :, :, :, :] = time_filter(i,self.grid.Nt) * self.alldata[i, :, :, :, :]

    def init_values(self):
        self.analize_field(self.filenames[0])

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
        if self.grid.Nt > 1:
            self.grid.dt = self.grid.t[1] - self.grid.t[0]
        else:
            self.grid.dt = 0
        self.grid.Lx = self.grid.dx * self.grid.Nx
        self.grid.Ly = self.grid.dy * self.grid.Ny
        self.grid.Lz = self.grid.dz * self.grid.Nz
        self.grid.Lt = self.grid.dt * self.grid.Nt
        self.grid.totNr = self.grid.Nt * self.grid.Nx * self.grid.Ny * self.grid.Nz

    def print_parameters(self):
        print ("Start Analysis:")
        print ("base name = %s" %self.basename)
        print ("SIZE: [ Lx, Ly, Lz ] = [ %.3f, %.3f, %.3f ]" % (self.grid.Lx, self.grid.Ly, self.grid.Lz))
        print (" Np : [ Nx, Ny, Nz ] = [ %d, %d, %d ]" % (self.grid.Nx, self.grid.Ny, self.grid.Nz))
        print ("time span = [%.1f:%.1f]     Nt = %d" % (self.tmin, self.tmax, self.grid.Nt))

    def set_frequency(self):

        self.grid.kx = np.fft.fftfreq(self.grid.Nx,d=self.grid.dx)
        self.grid.ky = np.fft.fftfreq(self.grid.Ny,d=self.grid.dy)
        self.grid.kt = np.fft.fftfreq(self.grid.Nt,d=self.grid.dt)

        self.grid.kx = np.fft.fftshift(self.grid.kx)
        self.grid.ky = np.fft.fftshift(self.grid.ky)
        #self.grid.kt = np.fft.fftshift(self.grid.kt)

        self.grid.Nkx = self.grid.kx.size
        self.grid.Nky = self.grid.ky.size
        self.grid.Nkt = self.grid.kt.size
        self.grid.totNk = self.grid.Nkx * self.grid.Nky * self.grid.kt.size

        self.grid.dkx = 1 / self.grid.Lx
        self.grid.dky = 1 / self.grid.Ly
        self.grid.dkt = 1 / self.grid.Lt

        self.grid.Lkx = self.grid.dkx * (self.grid.Nx / 2)
        self.grid.Lky = self.grid.dky * (self.grid.Ny / 2)
        self.grid.Lkt = self.grid.dkt * (self.grid.Nt / 2)

        self.grid.kz = 0
        self.grid.dkz = 0
        self.grid.Lkz = self.grid.dkz * (self.grid.Nz / 2)


    def print_frequency(self):
        print ("SIZE: [ Lkx, Lky, Lkz ] = [ %6.3f, %6.3f, %6.3f ]" % (self.grid.Lkx, self.grid.Lky, self.grid.Lkz))
        print (" dk : [ dkx, dky, dkz ] = [ %6.3f, %6.3f, %6.3f ]" % (self.grid.dkx, self.grid.dky, self.grid.dkz))
        print ("max_freq = %.2f     domega = %.2f" % (self.grid.Lkt, self.grid.dkt))


    def printEnergy(self, zposition, comp):
        factor = 8.0/(3.0*self.grid.Nt)*(self.grid.dx * self.grid.dy*2) / (8*math.pi)

        dataselect = self.alldata[:, zposition, :, :, comp]
        self.totalEnergyFunction = np.tensordot(dataselect,dataselect, axes=3)
        self.totalEnergyFunction = factor * self.totalEnergyFunction
        self.totalEnergyFFT = np.tensordot(self.trasf3D, self.trasf3D.conjugate(),axes=3)
        self.totalEnergyFFT = factor * np.absolute(self.totalEnergyFFT)
        print "totalEnergyFunction=" + str(self.totalEnergyFunction) + "   totalEnergyFFT=" + str(self.totalEnergyFFT)
        print "Efft/Einit=" + str(self.totalEnergyFFT/self.totalEnergyFunction)


    def do_fft(self, zposition, comp):
        print ("ready for fft3D...")
        dataselect = self.alldata[:, zposition, :, :, comp]
        self.trasf3D = np.fft.ifftn(dataselect, axes=(0,))
        self.trasf3D = self.trasf3D*math.sqrt(self.grid.Nt)
        self.trasf3D = np.fft.fftn(self.trasf3D, axes=(1,2))
        self.trasf3D = self.trasf3D/math.sqrt(self.grid.Nx*self.grid.Ny)

        self.shiftedTrasf3D = np.fft.fftshift(self.trasf3D, axes=(1,2))

        print ("DONE fft3D")

    def do_inversefft(self):
        print ("ready for ifft3D...")

        self.mynewData = np.fft.ifftn(self.trasf3D, axes=(1,2))
        self.mynewData = self.mynewData*math.sqrt(self.grid.Nx*self.grid.Ny)
        self.mynewData = np.fft.fftn(self.mynewData, axes=(0,))
        self.mynewData = self.mynewData/math.sqrt(self.grid.Nt)
        print ("DONE ifft3D")

    def save3Dfft(self):
        name = ("%s-3D-fft.txt" % (self.basename,))

        f1 = open(name, 'w')
        for t in range(0, self.grid.Nkt):
            for j in range(0, self.grid.Nky):
                for i in range(0, self.grid.Nkx):
                    f1.write("%e, %e, %e, %e\n" % (self.grid.kx[i], self.grid.ky[j], self.grid.kt[t], np.absolute(self.shiftedTrasf3D[t,j, i]) ) )
        # np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
        f1.close()

    def saveVTK3Dfft(self, varname, appendix=""):
        name = ("%s-fft%s.vtk" % (self.basename,appendix))
        factor = 2.*np.pi*0.00318
        writeLengthNt = (self.grid.Nkt/2)
        totPts = self.grid.Nkx * self.grid.Nky * writeLengthNt
        f1 = open(name, 'w')
        f1.write("# vtk DataFile Version 2.0\n")
        f1.write("titolo mio\n")
        f1.write("ASCII\n")
        f1.write("DATASET STRUCTURED_POINTS\n")
        f1.write("DIMENSIONS %d %d %d\n" % (self.grid.Nkx, self.grid.Nky, writeLengthNt))
        f1.write("ORIGIN %f %f %f\n" %(self.grid.kx[0]*factor, self.grid.ky[0]*factor, self.grid.kt[0]))
        f1.write("SPACING %f %f %f\n" % (self.grid.dkx*factor, self.grid.dky*factor, self.grid.dkt))
        f1.write("POINT_DATA %d\n" %totPts)
        f1.write("SCALARS %s float\n" %(varname))
        f1.write("LOOKUP_TABLE default\n")

        for t in range(0, writeLengthNt):
            for j in range(0, self.grid.Nky):
                for i in range(0, self.grid.Nkx):
                    f1.write("%f\n" %  np.absolute(self.shiftedTrasf3D[t, j, i]))
        # np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
        f1.close()


    def saveKxOmega(self):
        name = ("%s-kx-omega.txt" % (self.basename,))
        writeLengthNt = (self.grid.Nkt/2)

        f1 = open(name, 'w')
        for t in range(0, writeLengthNt):
            for i in range(0, self.grid.Nkx):
                f1.write("%e  %e %e\n" % (self.grid.kx[i], self.grid.kt[t], np.absolute(self.shiftedTrasf3D[t,(self.grid.ky.size)/2, i]) ) )
            f1.write("\n")
        # np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
        f1.close()

    def saveKxKyatomega(self, omegain):

        ifreqin = int(round(omegain/self.grid.dkt))
        for deltaifreq in range (-2,3,1):
            ifreqout=ifreqin+deltaifreq
            if(ifreqout>=0):
                omegaout = self.grid.kt[ifreqout]
                name = ("%s-kx-ky-at-omega%4.2f.txt" % (self.basename, omegaout))
                f1 = open(name, 'w')
                for j in range(0, self.grid.ky.size):
                    for i in range(0, self.grid.Nkx):
                        f1.write("%e, %e, %e\n" % (self.grid.kx[i], self.grid.ky[j], np.absolute(self.shiftedTrasf3D[ifreqout,j, i]) ) )
                    f1.write("\n")
                f1.close()

    def saveNewData(self, timein):

        itime = int(round(timein/self.grid.dt))
        name = ("%s-x-y.txt" % (self.basename))
        f1 = open(name, 'w')
        for j in range(0, self.grid.Ny):
            for i in range(0, self.grid.Nx):
                f1.write("%e, %e, %e\n" % (self.x[i], self.y[j], (self.mynewData[itime,j, i]) ) )
            f1.write("\n")
        f1.close()

    def setAllToZeroExceptOmega(self, omegain):

        ifreqin = int(round(omegain/self.grid.dkt))
        print "ifreqin=" + str(ifreqin) + " omegain=" + str(omegain)
        omegaout = self.grid.kt[ifreqin]
        print " omegaout=" + str(omegaout)
        for t in range(0, self.grid.Nkt/2):
            if abs(t-ifreqin)>1:
                self.trasf3D[t,:, :] = 0
                if t > 0:
                    self.trasf3D[self.grid.Nkt-t,:, :] = 0
        self.shiftedTrasf3D = np.fft.fftshift(self.trasf3D, axes=(1,2))

    def setToZeroOmega(self, omegain):

        ifreqin = int(round(omegain/self.grid.dkt))
        for t in range(0, self.grid.Nkt/2):
            if abs(t-ifreqin)<=1:
                self.shiftedTrasf3D[t,:, :] = 0
                self.trasf3D[t,:, :] = 0
                self.trasf3D[-t,:, :] = 0
        self.shiftedTrasf3D = np.fft.fftshift(self.trasf3D, axes=(1,2))



# def run():
#     print 'Number of arguments:', len(sys.argv), 'arguments.'
#     print 'Argument List:', str(sys.argv)
#     path = os.getcwd()
#     # f        = open(os.path.join(path,'E_FIELD_000.000.bin.000'),'rb')
#     component = 0
#     basetime = 0
#     deltaT = 15
#     sub = 5
#     endtime = basetime + deltaT
#     myname = "E_FIELD_"

#     for start in range(260, 350, 10):
#         print("analizzo il numero %d\n" % start)
#         basetime = start
#         endtime = basetime + deltaT
#         myanalysis = fieldAnalysis(basetime, endtime, substeps=sub, base=myname, end=".bin.000")
#         myanalysis.collect_data()
#         myanalysis.setWindow()
#         myanalysis.do_fft(zposition=0, comp=component)
#         myanalysis.saveVTK3Dfft()
#         myanalysis.saveKxOmega()
#run()

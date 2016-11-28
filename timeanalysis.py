###!/usr/bin/python
##loading shell commands
import os
import os.path
import glob
import sys
import shutil
import time
from scipy.fftpack import fftn
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
    return 0.5*(1-math.cos(2*np.pi*t/tmax))


class FieldAnalysis:
    """ la mia classe"""

    def __init__(self, tmin, tmax, substeps, base, end, reflect):
        self.base = base
        self.end = end
        self.tmin = int(tmin)
        self.tmax = int(tmax)
        self.substeps = int(substeps)
        self.reflect=reflect
        print "reflect:"
        print reflect
        if self.reflect:
            self.tmax=2*self.tmax
        if 1000 % substeps:
            print "ERROR: wrong number of timesteps!"
            sys.exit()
        self.Nt = substeps * (self.tmax - self.tmin)
        self.filenames = list()
        self.times = np.zeros((self.Nt))
        self.chose_file()

        self.init_values()
        self.set_frequency()
        self.basename = ("%s%03d" % (self.base, (self.tmin+self.tmin)/2 ))

    def analize_field(self, filename):
        f = open(filename, 'rb')
        # - endianness 0:small - 1:big -#

        endianness = struct.unpack('i', f.read(4))[0]

        # - global grid dim -#
        nx = struct.unpack('i', f.read(4))[0]
        ny = struct.unpack('i', f.read(4))[0]
        nz = struct.unpack('i', f.read(4))[0]
        self.Nx = nx
        self.Ny = ny
        self.Nz = nz
        ntot = nx * ny * nz
        # - processor grid -#
        npx = struct.unpack('i', f.read(4))[0]
        npy = struct.unpack('i', f.read(4))[0]
        npz = struct.unpack('i', f.read(4))[0]

        nproc = npx * npy * npz

        # - field components -#
        nc = struct.unpack('i', f.read(4))[0]
        self.Nx = nx
        self.Ny = ny
        self.Nz = nz
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

    # def collect_file(self):
    #     name = "E_FIELD_000.000.bin.000"
    #     for i in os.listdir(os.getcwd()):
    #         if i.endswith(".bin.000") and i.startswith("E_FIELD_"):
    #             # print i
    #             self.Nt += 1
    #             continue
    #         else:
    #             continue

    def chose_file(self):
        base = self.base
        end = self.end
        tindex = 0
        if self.reflect:
            for i in range(self.tmin, (self.tmax)/2):
                print i
                for v in range(0, self.substeps):
                    t = i + v * 1.0 / self.substeps
                    self.times[tindex] = t
                    name = base + ("%07.3f" % t) + end
                    self.filenames.append(name)
                    print name
                    tindex += 1
            print "PAUSA"
            for i in range((self.tmax)/2, self.tmin, -1):
                print i
                for v in range(0, self.substeps):
                    t = i + v * 1.0 / self.substeps
                    self.times[tindex] = t
                    t = i - v * 1.0 / self.substeps
                    name = base + ("%07.3f" % t) + end
                    self.filenames.append(name)
                    print name
                    tindex += 1
        else:
            for i in range(self.tmin, self.tmax):
                for v in range(0, self.substeps):
                    t = i + v * 1.0 / self.substeps

                    self.times[tindex] = t
                    name = base + ("%03d.%03d" % (i, v * 1000 / self.substeps)) + end
                    self.filenames.append(name)
                    # print name
                    tindex += 1
        if (self.Nt != len(self.filenames)):
            print ("error self.Nt=%d len(self.filenames))=%d" %(self.Nt, len(self.filenames)))

    def collect_data(self):
        self.alldata = np.zeros((self.Nt, self.Nz, self.Ny, self.Nx, self.Nc))
        for i in range(0, self.Nt):
            # print self.filenames[i]
            F = self.analize_field(self.filenames[i])
            self.alldata[i, :, :, :, :] = F[:, :, :, :]

    def setWindow(self):
        for i in range(0, self.Nt):
            self.alldata[i, :, :, :, :] = time_filter(i,self.Nt) * self.alldata[i, :, :, :, :]

    def init_values(self):
        self.analize_field(self.filenames[0])

        if self.Nx > 1:
            self.dx = self.x[1] - self.x[0]
        else:
            self.dx = 0
        if self.Ny > 1:
            self.dy = self.y[1] - self.y[0]
        else:
            self.dy = 0
        if self.Nz > 1:
            self.dz = self.z[1] - self.z[0]
        else:
            self.dz = 0
        if self.Nt > 1:
            self.dt = self.times[1] - self.times[0]
        else:
            self.dt = 0
        self.Lx = self.dx * self.Nx
        self.Ly = self.dy * self.Ny
        self.Lz = self.dz * self.Nz
        self.Lt = self.dt * self.Nt
        self.print_parameters()

    def print_parameters(self):
        print ("SIZE: [ Lx, Ly, Lz ] = [ %f, %f, %f ]" % (self.Lx, self.Ly, self.Lz))
        print (" Np : [ Nx, Ny, Nz ] = [ %d, %d, %d ]" % (self.Nx, self.Ny, self.Nz))
        print ("time span = [%f:%f]     Nt = %d" % (self.tmin, self.tmax, self.Nt))

    def set_frequency(self):
        
        self.kx = np.fft.fftfreq(self.Nx,d=self.dx)
        self.ky = np.fft.fftfreq(self.Ny,d=self.dy)
        kt = np.fft.fftfreq(self.Nt,d=self.dt)

        self.kx = np.fft.fftshift(self.kx)
        self.ky = np.fft.fftshift(self.ky)
        kt = np.fft.fftshift(kt)

        self.kt=np.zeros(kt.size/2)
        self.kt=-kt[kt.size/2:1:-1]
        
        self.dkx = 1 / self.Lx
        self.dky = 1 / self.Ly
        self.dkt = 1 / self.Lt
        
        self.Lkx = self.dkx * (self.Nx / 2)
        self.Lky = self.dky * (self.Ny / 2)
        self.Lkt = self.dkt * (self.Nt / 2)
        
        self.kz = 0
        self.dkz = 0
        self.Lkz = self.dkz * (self.Nz / 2)
        
        self.print_frequency()

        

    def print_frequency(self):
        print ("SIZE: [ Lkx, Lky, Lkz ] = [ %f, %f, %f ]" % (self.Lkx, self.Lky, self.Lkz))
        print (" dk : [ dkx, dky, dkz ] = [ %f, %f, %f ]" % (self.dkx, self.dky, self.dkz))
        print ("max_freq = %f     domega = %f" % (self.Lkt, self.dkt))

    def do_fft(self, zposition, comp):
        #trasf3D = np.zeros((self.Nt, self.Ny, self.Nx))
        #self.mytrasf3D = np.zeros((self.Nt/2, self.Ny, self.Nx))
        print ("ready for fft3D...")
        trasf3D = np.fft.fftn(self.alldata[:, zposition, :, :, comp])
        trasf3D = np.fft.fftshift(trasf3D)
        self.mytrasf3D =  trasf3D[self.Nt/2:1:-1,:,:]
        print ("DONE fft3D")
   
    def save3Dfft(self):
        name = ("%s-3D-fft.txt" % (self.basename,))

        f1 = open(name, 'w')
        for t in range(0, self.kt.size):
            for j in range(0, self.ky.size):
                for i in range(0, self.kx.size):
                    f1.write("%e, %e, %e, %e\n" % (self.kx[i], self.ky[j], self.kt[t], np.absolute(self.mytrasf3D[t,j, i]) ) )
        # np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
        f1.close()

    def saveVTK3Dfft(self):
        name = ("%s--fft.vtk" % (self.basename,))
        factor = 2.*np.pi*0.00318
        totPts = self.kx.size * self.ky.size * self.kt.size
        f1 = open(name, 'w')
        f1.write("# vtk DataFile Version 2.0\n")
        f1.write("titolo mio\n")
        f1.write("ASCII\n")
        f1.write("DATASET STRUCTURED_POINTS\n")
        f1.write("DIMENSIONS %d %d %d\n" % (self.kx.size, self.ky.size, self.kt.size))
        f1.write("ORIGIN %f %f %f\n" %(self.kx[0]*factor, self.ky[0]*factor, self.kt[0]))
        f1.write("SPACING %f %f %f\n" % (self.dkx*factor, self.dky*factor, self.dkt))
        f1.write("POINT_DATA %d\n" %totPts)
        f1.write("SCALARS fftBz float\n")
        f1.write("LOOKUP_TABLE default\n")
        
        
        for t in range(0, self.kt.size):
            for j in range(0, self.ky.size):
                for i in range(0, self.kx.size):
                    f1.write("%f\n" %  np.absolute(self.mytrasf3D[t, j, i]))
        # np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
        f1.close()


    def saveKxOmega(self):
        name = ("%s-kx-omega.txt" % (self.basename,))

        f1 = open(name, 'w')
        for t in range(0, self.kt.size):
            for i in range(0, self.kx.size):
                f1.write("%e, %e, %e\n" % (self.kx[i], self.kt[t], np.absolute(self.mytrasf3D[t,(self.ky.size)/2, i]) ) )
        # np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
        f1.close()


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


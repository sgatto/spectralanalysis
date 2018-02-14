########################################################################################################################
#                                                                                                                      #
#                            SPATIAL / TEMPORAL & SPATIO-TEMPORAL FOURIER TRANSFORM TOOLS                              #
#                                                                                                                      #
########################################################################################################################
import numpy as np

def to_w(t):
    w = np.fft.rfftfreq(t.size,d=(t[1]-t[0]) / (2.*np.pi))
    return w

def to_t(w, out_shape=None):
    nt   = 2*(w.size-1)+(out_shape[0]%2 if out_shape else 0 )
    t  = np.fft.fftfreq(nt,d=(w[1]-w[0])/(2.*np.pi))
    t  = np.fft.fftshift(t)
    t -= t[0]
    return t

def to_k(x):
    k = np.fft.fftfreq(x.size,d=(x[1]-x[0]) / (2.*np.pi))
    k = np.fft.fftshift(k)
    return k

def to_minus_k(x):
    k = -to_k(x)[::-1]
    return k

def to_x(k):
    x  = np.fft.fftfreq(k.size,d=(k[1]-k[0])/(2.*np.pi))
    x  = np.fft.fftshift(x)
    x -= x[0]
    return x


# TEMPORAL FOURIER TRANSFORM -------------------------------------------------------------------------------------------
def dFFT_t(t,f):
    w = to_w(t)
    g = np.fft.rfft(f,norm="ortho")
    return w, g

def iFFT_t(w,g, out_shape=None):
    t = to_t(w, out_shape )
    f = np.fft.irfft(g,norm="ortho", s=out_shape)
    return t, f

# 1D SPATIAL FOURIER TRANSFORM -----------------------------------------------------------------------------------------
def dFFT_x1(x,f):
    k = to_k(x)  
    g = np.fft.fftshift(np.fft.fft(f,norm="ortho"))
    return k, g

def iFFT_x1(k,g):
    x = to_x(k)
    f = np.fft.ifft(np.fft.ifftshift(g),norm="ortho")
    return x, f


# 2D SPATIAL FOURIER TRANSFORM -----------------------------------------------------------------------------------------
def dFFT_x2(x,y,f):
    kx = to_k(x)
    ky = to_k(y)
    g  = np.fft.fftshift( np.fft.fft2(f, axes=(1,0), norm="ortho") )
    return kx, ky, g

def iFFT_x2(kx,ky,g):
    x = to_x(kx)
    y = to_x(ky)
    f = np.fft.ifft2( np.fft.ifftshift(g), axes=(1,0), norm="ortho" )
    return x, y, f


# 3D SPATIAL FOURIER TRANSFORM -----------------------------------------------------------------------------------------
def dFFT_x3(x,y,z,f):
    kx = to_k(x)
    ky = to_k(y)
    kz = to_k(z)
    g  = np.fft.fftshift( np.fft.fftn(f, axes=(2,1,0), norm="ortho") )
    return kx, ky, kz, g

def iFFT_x3(kx,ky,kz,g):
    f = np.fft.ifftn( np.fft.ifftshift(g), axes=(2,1,0), norm="ortho")
    x = to_x(kx)
    y = to_x(ky)
    z = to_x(kz)
    return x, y, z, f


# 1D SPATIO-TEMPORAL FOURIER TRANSFORM ---------------------------------------------------------------------------------
def dFFT_tx1(t,x,f):
    w  = to_w(t)
    kx = to_minus_k(x)
    g  = np.fft.fftshift( np.fft.rfft2(f, norm="ortho", axes=(1,0)) , axes=(1) )[:,::-1]
    return w, kx, g

def iFFT_tx1(w,kx,g, out_shape=None):
    t = to_t(w, out_shape )
    x = to_x(kx)
    f = np.fft.irfft2( np.fft.ifftshift(g[:,::-1] , axes=(1) ), norm="ortho", axes=(1,0), s=out_shape)
    return t, x, f


# 2D SPATIO-TEMPORAL FOURIER TRANSFORM ---------------------------------------------------------------------------------
def dFFT_tx2(t,x,y,f):
    w  = to_w(t)
    kx = to_minus_k(x)
    ky = to_minus_k(y)
    g  = np.fft.fftshift( np.fft.rfftn(f, norm="ortho", axes=(2,1,0)) , axes=(1,2) )[:,::-1,::-1]
    return w, kx, ky, g

def iFFT_tx2(w,kx,ky,g, out_shape=None):
    t = to_t(w, out_shape )
    x = to_x(kx)
    y = to_x(ky)
    f = np.fft.irfftn( np.fft.ifftshift(g[:,::-1,::-1] , axes=(1,2) ), norm="ortho", axes=(2,1,0), s=out_shape)
    return t, x, y, f

# 3D SPATIO-TEMPORAL FOURIER TRANSFORM ---------------------------------------------------------------------------------

def dFFT_tx3(t,x,y,z,f):
    w  = to_w(t)
    kx = to_minus_k(x)
    ky = to_minus_k(y)
    kz = to_minus_k(z)
    g  = np.fft.fftshift( np.fft.rfftn(f, norm="ortho", axes=(3,2,1,0)) , axes=(1,2,3) )[:,::-1,::-1,::-1]
    return w, kx, ky, kz, g

def iFFT_tx3(w,kx,ky,kz,g, out_shape=None):
    t = to_t(w, out_shape )
    x    = to_x(kx)
    y    = to_x(ky)
    z    = to_x(kz)
    f = np.fft.irfftn( np.fft.ifftshift(g[:,::-1,::-1,::-1] , axes=(1,2,3) ), norm="ortho", axes=(3,2,1,0), s=out_shape)
    return t, x, y, z, f


# def scrivi_file2D(x, y, f, name):
#     f1 = open(name, 'w')
#     for i in range(0, x.size):
#         for j in range(0, y.size):
#             f1.write("%e, %e, %e\n" % (y[j], x[i], np.absolute(f[i,j])) )
#         f1.write("\n")
#     f1.close()
# 
# nt, nx, ny, ft, fx, fy = 256,256,256, 6, 10, 12
# t = np.arange(nt)*2*np.pi/nt
# x = np.arange(nx)*2*np.pi/nx
# y = np.arange(ny)*2*np.pi/ny
# A = np.fromfunction(lambda i,j,k: np.sin(2*np.pi*(ft*i/nt-(fx*j/nx+fy*k/ny))) , (nt,nx,ny), dtype=float)
# print A.shape
# 
# w, kx, ky, g = dFFT_tx2(t,x,y,A)
# print w.shape
# print kx.shape
# print ky.shape
# print g.shape
# print np.absolute(g) > 0.01 
# 
# t2, x2, y2, A2 = iFFT_tx2(w, kx, ky, g)
# 
# print "t2",t2.shape
# print "x2",x2.shape
# print "y2",y2.shape
# print "A2",A2.shape
# scrivi_file2D(w, k, g[:,:,140] , "pippo12.txt")



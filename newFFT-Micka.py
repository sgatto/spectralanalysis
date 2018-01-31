########################################################################################################################
#                                                                                                                      #
#                                                                                                                      #
#                            SPATIAL / TEMPORAL & SPATIO-TEMPORAL FOURIER TRANSFORM TOOLS                              #
#                                                                                                                      #
#                                                                                                                      #
########################################################################################################################
import numpy as np


# CONSTANT DEFINITION --------------------------------------------------------------------------------------------------
_twopi    = 2.*np.pi
_sq_twopi = np.sqrt(_twopi)



# TEMPORAL FOURIER TRANSFORM -------------------------------------------------------------------------------------------

# FORWARD
def tFFT(t,ft):
    # compute dft & axis
    fw  = np.fft.rfft( np.concatenate((ft[0::-1],ft[-1:0:-1]),axis=0),norm="ortho")
    w   = np.fft.rfftfreq(t.size,d=(t[1]-t[0])/ _twopi) 
    # return axes & Fourier transform
    return w, fw

# BACKWARD
def tFFT_inv(w,fw):
    # compute inverse dft & axis
    ft  = np.fft.irfft(fw,norm="ortho")
    t   = np.fft.rfftfreq(w.size,d=(w[1]-w[0])/ _twopi)
    ft = np.concatenate((ft[0::-1],ft[-1:0:-1]),axis=0)
    # return axes & Fourier transform
    return t, ft

# 1D SPATIAL FOURIER TRANSFORM -----------------------------------------------------------------------------------------
# FORWARD
def sFFT1(x,f_x):
    # compute dft & axis
    f_k  = np.fft.rfft(f_x,norm="ortho")
    k    = np.fft.rfftfreq(x.size,d=(x[1]-x[0])/ _twopi)
    return k, f_k


# BACKWARD
def sFFT1_inv(k,f_k):
    # compute inverse dft & axis
    f_x  = np.fft.irfft(f_k,norm="ortho")
    x    = np.fft.rfftfreq(k.size,d=(k[1]-k[0]) / _twopi)
    return x, f_x



# 2D SPATIAL FOURIER TRANSFORM -----------------------------------------------------------------------------------------

# FORWARD
def sFFT2(x,y,f):
    # compute dft & axis
    fkk   = np.fft.fft2(f)
    kx    = np.fft.fftfreq(x.size,d=x[1]-x[0]) * _twopi
    ky    = np.fft.fftfreq(y.size,d=y[1]-y[0]) * _twopi
    Kx,Ky = np.meshgrid(kx,ky)
    # add normalisation & phase-factor
    fkk  *= (x[1]-x[0])*(y[1]-y[0])/_twopi * np.exp(-complex(0,1)*Kx*x[0]-complex(0,1)*Ky*y[0])
    # shift ordering
    kx   = np.fft.fftshift(kx)
    ky   = np.fft.fftshift(ky)
    fkk  = np.fft.fftshift(fkk)
    return kx, ky, fkk

# BACKWARD
def sFFT2_inv(kx,ky,fkk):
    # compute inverse dft & axis
    f    = np.fft.ifft2(fkk)
    x    = np.fft.fftfreq(kx.size,d=kx[1]-kx[0]) * _twopi
    y    = np.fft.fftfreq(ky.size,d=ky[1]-ky[0]) * _twopi
    X,Y  = np.meshgrid(x,y)
    # add normalisation & phase-factor
    f   *= kx.size*ky.size*(kx[1]-kx[0])*(ky[1]-ky[0])/_twopi * np.exp(-complex(0,1)*X*kx[0]-complex(0,1)*Y*ky[0])
    # shift ordering
    x    = np.fft.fftshift(x)
    y    = np.fft.fftshift(y)
    f_xy = np.fft.fftshift(f)
    return x, y, f_xy


####################################################################################### UNDER DEVELOPMENT (NOT YET WORKING)

# 1D SPATIO-TEMPORAL FOURIER TRANSFORM ---------------------------------------------------------------------------------

# FORWARD
def STfft1(t,x,f_tx):
    # compute dft & axis
    f_wk  = np.fft.ifftn(f_tx, axes=(0,))
    w     = np.fft.fftfreq(t.size,d=t[1]-t[0]) * _twopi
    f_wk  = np.fft.fftn(f_wk, axes=(1,))
    k     = np.fft.fftfreq(x.size,d=x[1]-x[0]) * _twopi
    W,K   = np.meshgrid(w,k)
    # add normalisation & phase-factor
    f_wk *= (t[1]-t[0])*(x[1]-x[0])/_twopi # * np.exp(complex(0,1)*W*t[0]-complex(0,1)*K*x[0])
    # shift ordering
    #w    = np.fft.fftshift(w)
    k    = np.fft.fftshift(k)
    f_wk = np.fft.fftshift(f_wk, axes=(1,))
    # return axes & Fourier transform
    return w, k, f_wk

# FORWARD
def STfft2(t,x,y,f_tx):
    # compute dft & axis
    f_wk  = np.fft.rfftn(np.flip(f_tx,axis=0), axes=(1,2,0), norm="ortho")
    w     = np.fft.rfftfreq(t.size,d=t[1]-t[0]) * _twopi
    kx    = np.fft.fftfreq(x.size,d=x[1]-x[0]) * _twopi
    ky    = np.fft.fftfreq(y.size,d=y[1]-y[0]) * _twopi
    # shift ordering
    kx   = np.fft.fftshift(kx)
    ky   = np.fft.fftshift(ky)
    f_wk = np.fft.fftshift(f_wk, axes=(1,2))
    # return axes & Fourier transform
    return w, kx, ky, f_wk

# # BACKWARD
# def iTfft(w,f_w):
#     # compute inverse dft & axis
#     f_t  = np.fft.ifft(f_w)
#     t    = np.fft.fftfreq(w.size,d=w[0]-w[1]) * _twopi
#     # add normalisation & phase-factor
#     f_t *= w.size*(w[1]-w[0])/_sq_twopi * np.exp(-complex(0,1)*t*w[0])
#     # shift ordering
#     t   = np.fft.fftshift(t)
#     f_t = np.fft.fftshift(f_t)
#     t   = t[::-1]
#     f_t = f_t[::-1]
#     # return axes & Fourier transform
#     return t, f_t

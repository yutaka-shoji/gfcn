#!/usr/bin/env python3

import numpy as np
from scipy import pi
import ctypes
import time

class gfunc(object):
    def __init__(self, model='ils'):
        self.model = model

        # INFINITE LINE SOURCE MODEL
        if model == 'ils':
            self.lib = np.ctypeslib.load_library('libils.so', '.')
            self.lib.calc_responsefcn.argtypes = [
                    ctypes.POINTER(ctypes.c_int),
                    np.ctypeslib.ndpointer(dtype=np.float64),
                    np.ctypeslib.ndpointer(dtype=np.float64)
                    ]
            self.lib.calc_responsefcn.restype = ctypes.c_void_p

        # INFINITE CYLINDER SOURCE MODEL
        if model == 'ics':
            self.lib = np.ctypeslib.load_library('libics.so', '.')
            self.lib.set_params.argtypes = [ ctypes.POINTER(ctypes.c_double) ]
            self.lib.set_params.restype = ctypes.c_void_p
            self.lib.calc_responsefcn.argtypes = [
                    ctypes.POINTER(ctypes.c_int),
                    np.ctypeslib.ndpointer(dtype=np.float64),
                    np.ctypeslib.ndpointer(dtype=np.float64)
                    ]
            self.lib.calc_responsefcn.restype = ctypes.c_void_p

        # COMPOSITE-MEDIUM INFINITE LINE SOURCE MODEL
        if model == 'cmils':
            self.lib = np.ctypeslib.load_library('libcmils.so', '.')
            self.lib.set_params.argtypes = [
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double)
                    ]
            self.lib.set_params.restype = ctypes.c_void_p
            self.lib.calc_responsefcn.argtypes = [
                    ctypes.POINTER(ctypes.c_int),
                    np.ctypeslib.ndpointer(dtype=np.float64),
                    np.ctypeslib.ndpointer(dtype=np.float64)
                    ]
            self.lib.calc_responsefcn.restype = ctypes.c_void_p

def main():
    r_b = 0.060                                 # (m)
    r_po = 0.032 / 2                            # (m)
    D = 0.035                                   # (m)
    RA = (D + r_po) / r_b
    RB = (D - r_po) / r_b
    RD = D / r_b
    R = r_b / r_b
    ks = 2.50                                   # (W m-1 K-1)
    Cs = 3.0e+6                                 # (J m-3 K-1)
    alpha_s = ks / Cs                           # (m2 s-1)
    kb = 1.80
    Cb = 3.8e+6
    #kb = ks                                    # (W m-1 K-1)
    #Cb = Cs                                    # (J m-3 K-1)
    alpha_b = kb / Cb                           # (m2 s-1)
    alpha = np.sqrt(alpha_b / alpha_s)
    k = ks / kb

    tmin = 60.0                                 # (s)
    tmax = 60.0 * 60 * 24 * 365 * 10            # (s)
    t = np.logspace(np.log10( tmin ), np.log10( tmax + 1 ), num=100, endpoint=True)

    # ILS-------------------------------------------------
    start = time.time()
    Fo_ils = alpha_s * t / (r_b*r_b)
    #Fo_ils = alpha_s * t / (D*D)
    gfunc_ils = gfunc('ils')
    f_ils = np.zeros_like(Fo_ils)
    n = ctypes.byref(ctypes.c_int(Fo_ils.size))
    gfunc_ils.lib.calc_responsefcn(n, Fo_ils, f_ils)
    G_ils = 1 / (4 * pi * ks) * f_ils
    print(f'  ils: {time.time() - start}')
    # ILS-------------------------------------------------

    # ICS-------------------------------------------------
    start = time.time()
    Fo_ics = alpha_s * t / (r_b*r_b)
    gfunc_ics = gfunc('ics')
    gfunc_ics.lib.set_params( ctypes.byref(ctypes.c_double(R)) )
    f_ics = np.zeros_like(Fo_ics)
    n = ctypes.byref(ctypes.c_int(Fo_ics.size))
    gfunc_ics.lib.calc_responsefcn(n, Fo_ics, f_ics)
    G_ics = 1 / (pi*pi * ks) * f_ics
    print(f'  ics: {time.time() - start}')
    # CMILS-------------------------------------------------

    # CMILS-------------------------------------------------
    start = time.time()
    Fo_cmils = alpha_b * t / (r_b*r_b)
    gfunc_cmils = gfunc('cmils')
    gfunc_cmils.lib.set_params(
            ctypes.byref(ctypes.c_double(RA)),
            ctypes.byref(ctypes.c_double(RB)),
            ctypes.byref(ctypes.c_double(RD)),
            ctypes.byref(ctypes.c_double(alpha)),
            ctypes.byref(ctypes.c_double(k))
            )
    # full
    f_cmils = np.zeros_like(Fo_cmils)
    n = ctypes.byref(ctypes.c_int(Fo_cmils.size))
    gfunc_cmils.lib.calc_responsefcn(n, Fo_cmils, f_cmils)
    G_cmils = 1 / (2 * pi * kb) * f_cmils
    print(f'cmils: {time.time() - start}')
    # CMILS-------------------------------------------------

    # plot
    import matplotlib.pyplot as plt
    plt.rcParams['font.size'] = 10
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'stix'
    fig = plt.figure(figsize=(5,5), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    sc = ax.scatter(
            t, G_ils,
            facecolor = 'none',
            edgecolor = 'black',
            marker = 'o',
            s = 20,
            lw = 0.5,
            label = 'ILS'
            )
    sc = ax.scatter(
            t, G_ics,
            color = 'black',
            marker = '+',
            s = 20,
            lw = 0.5,
            label = 'ICS'
            )
    sc = ax.scatter(
            t, G_cmils,
            facecolor = 'none',
            edgecolor = 'black',
            marker = '^',
            s = 20,
            lw = 0.5,
            label = 'CM-ILS (at U-tube wall)'
            )
    ax.set_xlabel('time')
    ax.set_ylabel('$G(t)$')
    ax.set_xscale('log')
    ax.set_title(
            r'$k_\mathrm{s} = %.2f$  $k_\mathrm{b} = %.2f$' % (ks, kb),
            fontsize=10
            )
    ax.set_xticks( [ 60, 60 * 60, 60 * 60 * 24, 60 * 60 * 24 * 365, 60 * 60 * 24 * 365 * 10 ] )
    ax.set_xticklabels( ['1 min', '1 hour', '1 day', '1 year', '10 year'] )
    ax.minorticks_off()
    plt.legend(frameon=False)
    plt.show()

if __name__ == '__main__':
    main()

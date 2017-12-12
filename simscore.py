# -*- coding: utf-8 -*-
"""
Created on Fri Dec 01 15:23:42 2017

Original from Yim et al. 2015.
Modified to use the .npz file system instead of the .txt from Yim.
"""

import numpy as np
import pylab

def simscore(file1, file2, delta, bin = 1.0, dur=100.0, ncell=2000):
    """Calculate the similarity score as in Yim et al. 2015"""
    d1 = np.loadtxt(file1)
    d2 = np.loadtxt(file2)
    x = np.zeros(int(ncell*dur/bin))
    y = np.zeros(int(ncell*dur/bin))

    for j in range(ncell):
        if d1.size == 2:
            s1 = np.array(d1[0]*(d1[1]==j))
        else:
            s1 = d1[d1[:,1]==j,0]
        if d2.size == 2:
            s2 = np.array(d2[0]*(d2[1]==j))
        else:
            s2 = d2[d2[:,1]==j,0]

        kern = np.append(np.arange(delta/bin),np.arange(delta/bin,-1,-1))
        ts1,dump = pylab.histogram(s1,np.arange(0.,dur+bin,bin))
        ts2,dump = pylab.histogram(s2,np.arange(0.,dur+bin,bin))
        x[j*dur/bin:(j+1)*dur/bin] = np.convolve(ts1,kern,'same')
        y[j*dur/bin:(j+1)*dur/bin] = np.convolve(ts2,kern,'same')

    x = x - pylab.mean(x)
    y = y - pylab.mean(y)
    cor = sum(x*y)/(len(x)*pylab.std(x)*pylab.std(y))

    return cor
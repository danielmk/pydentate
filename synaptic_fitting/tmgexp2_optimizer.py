#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# PEP8
"""
Created on Thu Nov  7 14:48:40 2019

@author: can
"""
from scipy.optimize import minimize
import numpy as np
import time
from tmgexp2_simulator import simulate
from peak_finder_tmgexp2 import peakfinder
begin = time.time()


# Data path msut be given here either for gcmc or gcin, or another dataset
# unused one was commented out
data_path_gcmc = "/home/can/Downloads/gc_to_mc"
#data_path_gcin = "/home/can/Downloads/gc_to_in"
peaks = peakfinder(data_path_gcmc)[0]
#peaks = peakfinder(data_path_gcin)[0]
taus = []
'LOSS function'
# with Mean Square Error for each freq and then avarage of all'


def loss(x):
    tau_facil, tau_rec, u0 = x
    # taus and u0 are collected for each iteration of loss func
    taus.append(x)
    output1hz = simulate(x[0], x[1], 1, x[2])[0]
    Hz1 = peaks[0]
    output10hz = simulate(x[0], x[1], 10, x[2])[0]
    Hz10 = peaks[1]
    output30hz = simulate(x[0], x[1], 30, x[2])[0]
    Hz30 = peaks[2]
    output50hz = simulate(x[0], x[1], 50, x[2])[0]
    Hz50 = peaks[3]

    mse1 = (np.square(output1hz - Hz1)).mean(axis=None)
    mse10 = (np.square(output10hz - Hz10)).mean(axis=None)
    mse30 = (np.square(output30hz - Hz30)).mean(axis=None)
    mse50 = (np.square(output50hz - Hz50)).mean(axis=None)
    mse = (mse1 + mse10 + mse30 + mse50)/4

    return mse


'OPTIMIZATION'

# initial array to optimize is given here as x0
x0 = np.array([500, 0, 0.1])
# resulting MSE of loss function is minimized by Nelder-Mead Algorithm
res = minimize(loss, x0, method='Nelder-Mead')
# results are saved
np.savez('gc_to_mc_opt', loss = res)
# np.savez('gc_to_in_opt', loss = res)
# time that optimizaation takes 
end = time.time()
opt_time = end-begin
print (opt_time)

"""
RESULTS

res for granule cell to interneuron, order [tau_facil, tau_rec, u0]
 
Out[2]: 
 final_simplex: (array([[4.11189527e+03, 8.78187450e-03, 8.05591091e-02]
           fun: 0.6866318365387405
       message: 'Optimization terminated successfully.'
          nfev: 221
           nit: 111
        status: 0
       success: True
             x: array([4.11189527e+03, 8.78187450e-03, 8.05591091e-02])
             
             
             

res for granule cell to mossy cell, order [tau_facil, tau_rec, u0]

Out[2]: 
 final_simplex: (array([[3.82057408e+03, 3.26025349e+01, 7.78641198e-02]
           fun: 0.028299227264135973
       message: 'Optimization terminated successfully.'
          nfev: 284
           nit: 155
        status: 0
       success: True
             x: array([3.82057408e+03, 3.26025349e+01, 7.78641198e-02])

"""
"""
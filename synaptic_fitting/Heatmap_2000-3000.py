#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:36:56 2019

@author: barisckuru
"""

import numpy as np
import os
from tmgexp2_simulator import simulate
import matplotlib.pyplot as plt
import time

begin = time.time()

# PARAMETERS
freq_1 = []
freq_10 = []
freq_30 = []
freq_50 = []
load_1 = []
peaks = []
taus = []
norm = []

'LOAD THE EXPERIMENTAL DATA'

data_path = "/home/can/Downloads/gc_to_mc/"
'Open and load the data in the same directory for diff freqs'
for files in os.walk(data_path):
    for file in files[2]:
        if '1hz' in file:
            curr_path = data_path + file
            load_1 = np.load(curr_path)
            freq_1 = load_1['mean_arr']
        if '10hz' in file:
            curr_path = data_path + file
            load_1 = np.load(curr_path)
            freq_10 = load_1['mean_arr']
        if '30hz' in file:
            curr_path = data_path + file
            load_1 = np.load(curr_path)
            freq_30 = load_1['mean_arr']
        if '50hz' in file:
            curr_path = data_path + file
            load_1 = np.load(curr_path)
            freq_50 = load_1['mean_arr']
loads = [freq_1, freq_10, freq_30, freq_50]

'PEAK FINDER'
for i in range(len(loads)):
    data = loads[i]
    '''Data is current response, values are negative.
    For stimulus artifact, both positive and negative
    positive threshold was used to define the indices for stimuli'''
    indices = np.argwhere(data > 200)
    'Indices were shifted 40 dps, response without artifact in the beginning'
    indices = np.transpose(indices)[0] + 40
    'Data was reverted to compare it positive conductance values from sim'
    data = -data
    'One more indice was appended to create an interval for the last stimulus'
    indices = np.append(indices, indices[-1] + (indices[2] - indices[1]))
    'NORMALIZATION of data by the local max of first signal'
    # end is cut to eliminate stim artifact
    first_sig = data[indices[1]:indices[3]-60]
    first_peak = max(first_sig)
    data = data/first_peak
    data_cut = data[(indices[1]-10000):]
    norm.append(data_cut)

    '''Indices for peak finder, 2 idcs for one stimulus, 1 picked
    Shifted and selected idcs are now starting points
    Stop points were also defined by adding the length of syn response
    so stim artifact at the end was eliminated'''
    start = indices[1::2]
    stop = start + len(first_sig)
    indices = np.concatenate((start, stop))
    indices = np.sort(indices)

    '''Data were splitted with respect to indices
    local max was found for each part'''
    split_data = np.split(data, indices)
    split_data = split_data[1:len(split_data):2]
    split_data = np.array(split_data)
    peaks_data = np.amax(split_data, axis=1)
    peaks.append(peaks_data)


'LOSS func w Mean Square Error for each freq and then avarage'


def loss(x):
    tau_facil, tau_rec = x
    taus.append(x)
    u0 = 7.78641198e-02
    sampling = 0.5
    output1hz = simulate(x[0], x[1], 1, u0, sampling)[0]
    Hz1 = peaks[0]
    output10hz = simulate(x[0], x[1], 10, u0, sampling)[0]
    Hz10 = peaks[1]
    output30hz = simulate(x[0], x[1], 30, u0, sampling)[0]
    Hz30 = peaks[2]
    output50hz = simulate(x[0], x[1], 50, u0, sampling)[0]
    Hz50 = peaks[3]

    mse1 = (np.square(output1hz - Hz1)).mean(axis=None)
    mse10 = (np.square(output10hz - Hz10)).mean(axis=None)
    mse30 = (np.square(output30hz - Hz30)).mean(axis=None)
    mse50 = (np.square(output50hz - Hz50)).mean(axis=None)
    mse = (mse1 + mse10 + mse30 + mse50)/4

    return mse


# with 1hz, 770 hours
# w/o 1hz 189 hours
    


Z = []
pars = []
mat = np.zeros((501,101))
tau_facil = np.arange(2000, 3001, 2)
tau_rec = np.arange(0, 201, 2)
for i in tau_facil:
    for j in tau_rec:
        x1 = np.array([i,j])
        pars.append(x1)
        curr_loss = loss(x1)
        idc_facil = int((int(i)-2000)/2)
        idc_rec = int(int(j)/2)
        mat[idc_facil, idc_rec] = curr_loss


np.savez('heatmap_2k-3k', loss=mat)

end = time.time()
print('time(seconds): ', end-begin)



loss_load1 = np.load('heatmap_2k-3k.npz')
loss_load2 = np.load('heatmap_3k-4k.npz')
loss_load3 = np.load('heatmap_4k-5k.npz')
loss1 = loss_load1['loss']
loss2 = loss_load2['loss']
loss3 = loss_load3['loss']


loss1 = np.delete(loss1, 500, axis = 0)
loss2 = np.delete(loss2, 500, axis = 0)
loss_all= np.concatenate((loss1,loss2,loss3))

np.savez('heatmap_gc_to_mc', loss=loss_all)

loss_all=np.load('heatmap_gc_to_mc.npz')
log_all = np.log(loss_all)
plt.imshow(log_all, aspect='auto', cmap=plt.cm.cividis, 
           interpolation='none', extent=[0,200,2000,5000])

plt.title('Heatmap GC to MC')
plt.xlabel('tau recovery')
plt.ylabel('tau facilitation')
plt.colorbar(label='log of loss values')

end = time.time()
print('time(seconds): ', end-begin)
"""
CS = plt.contourf(loss_all, cmap=plt.cm.cividis)
plt.colorbar()  
plt.show()

"""

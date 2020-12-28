#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:31:34 2019

@author: barisckuru
"""

import numpy as np
import matplotlib.pyplot as plt

'Heatmap GC to MC'

loss_gcmc = np.load('heatmap_gc_to_mc.npz')['loss']
log_gcmc = np.log(loss_gcmc)
plt.imshow(log_gcmc, aspect='auto', cmap=plt.cm.cividis, 
           interpolation='none', extent=[0,200,2000,5000])
plt.title('Heatmap GC to MC')
plt.xlabel('tau recovery')
plt.ylabel('tau facilitation')
plt.colorbar(label='log of loss values')

'Heatmap GC to IN'

loss_gcin = np.load('heatmap_gc_to_in.npz')['loss']
log_gcin = np.log(loss_gcin)
plt.figure()
plt.imshow(log_gcin, aspect='auto', cmap=plt.cm.cividis, 
           interpolation='none', extent=[0,200,2000,5000])
plt.title('Heatmap GC to IN')
plt.xlabel('tau recovery')
plt.ylabel('tau facilitation')
plt.colorbar(label='log of loss values')

"""
concatenation of 3 parts

for GCMC

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

for GCIN

loss_load1 = np.load('gcin_2k-3k.npz')
loss_load2 = np.load('gcin_3k-4k.npz')
loss_load3 = np.load('gcin_4k-5k.npz')
loss1 = loss_load1['loss']
loss2 = loss_load2['loss']
loss3 = loss_load3['loss']


loss1 = np.delete(loss1, 500, axis = 0)
loss2 = np.delete(loss2, 500, axis = 0)
loss_all= np.concatenate((loss1,loss2,loss3))
np.savez('heatmap_gc_to_in', loss=loss_all)
"""

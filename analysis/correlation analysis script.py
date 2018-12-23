# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 15:35:07 2018

@author: spell
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.stats import pearsonr
import os


def correlation_analysis(a,b):      #with 2 input arrays, returns array of row-wise correlations. array shape must not be 0
    alist = []
    n = 0
    x, y = a.shape
    for n in range (x):
        a2 = a[n,:]
        b2 = b[n,:]
        c = scipy.stats.pearsonr(a2,b2)
        d = c[0]
        alist.append(d)
    matrix = np.array(alist)
    return matrix

def correlation_average (a):            # with one input array returns the average of that array ignoring 'nan'
        # g = np.nan_to_num(a)
        h = np.nanmean(a)
        return h

def filenames_list(a): #returns a lits with all the filenames. a is the directory
    path = a
    folder = os.fsencode(path)
    filenames = []
    for file in os.listdir(folder):
        filename = os.fsdecode(file)
        if filename.endswith(".npz"): # whatever file types you're using.
            filenames.append((a + "\\" + filename))
    filenames.sort()
    return filenames

def directory_correlation(f):                   #takes file list from filenames_list function.
    alist = []
    m = 0
    for n in range(len(f)):             
        for m in range(len(f)):         
            my_data_a = np.load(f[n])
            my_data_b = np.load(f[m])
            a = my_data_a['arr_0']
            b = my_data_b['arr_0']
            c = correlation_average(correlation_analysis(a,b))
            alist.append(c)
    x = np.array(alist)
    y = np.reshape(x, ((len(f)),(len(f))))
    return y


a = "Z:\\pyDentate\\pyDentateData\\pattern_separation_data_local_input_revised\\seed10000\\input_patterns\\"
b = filenames_list(a)
c = directory_correlation(b)
print(c)

# or : directory_correlation(filenames_list("C:\\Users\\spell\\Desktop\\DanielDataShannanigans\\input_patterns_seed_1000"))

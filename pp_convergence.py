# -*- coding: utf-8 -*-
"""
Created on Tue May 01 11:01:42 2018

@author: daniel
"""

import numpy as np
import matplotlib.pyplot as plt

innervation_pattern_gc = np.array([np.random.choice(400,20, replace = False) for x in range(2000)])
innervation_pattern_gc = innervation_pattern_gc.swapaxes(0,1)

PP_to_GCs = []
for x in range(0,400):
    PP_to_GCs.append(np.argwhere(innervation_pattern_gc == x)[:,1])

PP_to_GCs = np.array(PP_to_GCs)

plt.figure()
hist_one = plt.hist(np.hstack(PP_to_GCs[0:24].flat), bins=2000)
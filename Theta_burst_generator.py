# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 15:37:20 2018

@author: DanielM
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(10000)

temporal_patterns_inter_burst = np.random.poisson(100, (400,20)).cumsum(axis=1)[:,10:20] - 1000
temporal_patterns_inter_burst = np.repeat(temporal_patterns_inter_burst[:,:,np.newaxis],5,axis=2)

temporal_patterns_intra_burst = np.random.poisson(10, (4000,10,5,))[3600:4000,:,:].cumsum(axis=1)

test = temporal_patterns_inter_burst + temporal_patterns_intra_burst
test2 = test.reshape(400,50)
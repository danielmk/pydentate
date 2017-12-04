# -*- coding: utf-8 -*-
"""
Created on Fri Dec 01 15:47:52 2017

@author: Holger
"""

import numpy as np

def get_poisson_pattern(rate = 10, npatterns = 1):
    return np.random.poisson(10, (npatterns, 3)).cumsum()
    
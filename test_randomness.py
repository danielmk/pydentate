# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:03:55 2017

@author: DanielM
"""

import numpy as np
import matplotlib.pyplot as plt

def cell_picker(target_pool, divergence):

    pre_pop_spat_norm = np.arange(500,dtype = float) / (500 - 1)
    post_pop_spat_norm = np.arange(15, dtype = float) / (15 - 1)

    #pre_pop_spat_norm = ((np.arange(500,dtype = float)-249)  / 250)
    #post_pop_spat_norm = ((np.arange(15,dtype = float)-7) / 7)

    pre_cell_target = []

    for curr_cell_spat in pre_pop_spat_norm:
        post_pop_dist = np.absolute(post_pop_spat_norm - curr_cell_spat)
        sorted_idc = np.argsort(post_pop_dist)
        random_int = np.random.choice(sorted_idc[:target_pool], divergence, replace = False)
        pre_cell_target.append(random_int)

    return np.array(pre_cell_target)

if __name__ == '__main__':
    results = []
    for x in range(1000):
        results.append(cell_picker(3,1))
    results = np.array(results)


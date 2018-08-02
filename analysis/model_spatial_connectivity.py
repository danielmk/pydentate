# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 10:25:41 2018

@author: daniel
"""
import shelve
import numpy as np
import matplotlib.pyplot as plt

def pos(rad):
    """
    (x,y) position of a point on a circle with axis origin at (0,0)
    and radius 1.
    x = cx + r * cos(rad) -> x = cos(rad)
    y = cy + r * sin(rad) -> y = sin(rad)
    
    Returns a list of tuples that give the point of each radian passed.
    """
    x_arr = list(np.cos(rad))
    y_arr = list(np.sin(rad))
    
    return [(x_arr[idx],y_arr[idx]) for idx in range(len(x_arr))]

def euclidian_dist(p1,p2):
    """ p1 and p2 must both be of len 2 where p1 = (x1,y1); p2 = (x2,y2)"""
    return math.sqrt((p1[0] - p2[0])**2.0 + (p1[1] - p2[1])**2.0)

#Home PC
#directory = "C:\\Users\\daniel\\repos\\pyDentate\paradigm_pattern-separation_saves_2018-03-11\\"
#Office PC
#directory = "Y:\\DanielM\\023_Dentate Gyrus Model\\paradigm_spatial-inhibition\\"
#Dropbox
directory = "C:\\Users\\Daniel\\repos\\pyDentate\\"

file_name = "net_tunedrevexpmoresyn.TunedNetwork_data_paradigm_local-pattern-separation_run_scale_seed_000_1000_10000.pydd"

data = shelve.open(directory + file_name)

# Get to BasketCell Connection
BC_to_GC_targets = data['net_tunedrevexpmoresyn.TunedNetwork']['populations'][0]['connections'][27]['BasketCellPopulation to GranuleCellPopulation']['pre_cell_targets']

pre_pop_rad = (np.arange(24,dtype=float) / 24.0) * (2*np.pi)
post_pop_rad = (BC_to_GC_targets / 2000.0) * (2*np.pi)
#post_pop_rad = (np.arange(2000, dtype=float) / 2000) * (2*np.pi)

pre_pop_pos = pos(pre_pop_rad)
post_pop_pos = pos(post_pop_rad)

pre_pop_pos_arr = np.array(pre_pop_pos)
post_pop_pos_arr = np.array(post_pop_pos)

BC_to_GC_arc = np.array(BC_to_GC_targets, dtype = np.float)

for row in range(BC_to_GC_arc.shape[0]):
    for col in range(BC_to_GC_arc.shape[1]):
        dist1 = abs(pre_pop_rad[row] - post_pop_rad[row,col])
        dist2 = abs(dist1 - (2*np.pi))
        if dist1 <= dist2:
            BC_to_GC_arc[row,col] = dist1
        elif dist1 > dist2:
            BC_to_GC_arc[row,col] = dist2

#BC_to_GC_arc = np.absolute(BC_to_GC_arc)

plt.figure()
myhist = plt.hist((BC_to_GC_arc.flatten() / (2*np.pi))*2000, np.arange(-0.5,2000,1))
plt.xlabel("relative #GC")
plt.ylabel("# Synapses")
plt.title("net_globalrev_run_scale_seed_000_1000_10000 BC->GC connection")

plt.figure()
cum_hist = myhist[0].cumsum() / myhist[0].cumsum()[-1]
plt.plot(cum_hist)
plt.xlabel("relative #GC")
plt.ylabel("Cumulated fraction of synapses")
plt.title("net_globalrev_run_scale_seed_000_1000_10000 BC->GC connection")

def pos(rad):
    """
    (x,y) position of a point on a circle with axis origin at (0,0)
    and radius 1.
    x = cx + r * cos(rad) -> x = cos(rad)
    y = cy + r * sin(rad) -> y = sin(rad)
    
    Returns a list of tuples that give the point of each radian passed.
    """
    x_arr = list(np.cos(rad))
    y_arr = list(np.sin(rad))
    
    return [(x_arr[idx],y_arr[idx]) for idx in range(len(x_arr))]

def euclidian_dist(p1,p2):
    """ p1 and p2 must both be of len 2 where p1 = (x1,y1); p2 = (x2,y2)"""
    return math.sqrt((p1[0] - p2[0])**2.0 + (p1[1] - p2[1])**2.0)


"""
plt.figure()
plt.hist(BC_to_GC_targets, bins = 2000)
plt.xlabel("# GCs")
plt.ylabel("# incoming BC Synapses")
"""
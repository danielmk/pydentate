# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 15:37:20 2018

@author: DanielM
"""

import numpy as np
import matplotlib.pyplot as plt
from elephant import spike_train_generation, statistics
import quantities as pq

def poisson_generator(interval =100, t_stop = 1000.0, nr_spikes = 10,
                            nr_patterns = 400, numpy_seed = 10000):
    np.random.seed(numpy_seed)
    temporal_patterns_inter_burst = []
    for x in range(nr_patterns):
        spike_train = np.random.exponential(interval, (nr_spikes)).cumsum()
        temporal_patterns_inter_burst.append(spike_train[spike_train <= t_stop])

    temporal_patterns_inter_burst = np.array(temporal_patterns_inter_burst)

    return temporal_patterns_inter_burst

def poisson_burst_generator(interval =100, t_stop = 1000.0, nr_spikes = 10,
                            nr_patterns = 400, numpy_seed = 10000):
    

#temporal_patterns_inter_burst = np.repeat(temporal_patterns_inter_burst[:,:,np.newaxis],5,axis=2)
"""
temporal_patterns_intra_burst = np.random.poisson(10, (4000,10,5,))[3600:4000,:,:].cumsum(axis=1)

test = temporal_patterns_inter_burst + temporal_patterns_intra_burst
test2 = test.reshape(400,50)"""

if __name__ == '__main__':
    temporal_patterns_inter_burst = poisson_generator(numpy_seed = 10000)
    plt.figure()
    plt.eventplot(temporal_patterns_inter_burst)
    plt.title("numpy generated eventplot")

    plt.figure()
    plt.hist(np.hstack(temporal_patterns_inter_burst.flat), bins = int(np.hstack(temporal_patterns_inter_burst.flat).max()))
    plt.title("numpy generated time histogram")

    isis = []
    for x in temporal_patterns_inter_burst:
        isi = statistics.isi(x)
        isis.append(isi)

    isis = np.array(isis)
    
    plt.figure()
    plt.hist(np.hstack(isis.flat), bins = range(0,int(np.hstack(isis.flat).max()-1)))
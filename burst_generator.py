# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 15:37:20 2018

@author: DanielM
"""

import numpy as np
import matplotlib.pyplot as plt

def poisson_generator(interval =100,
                      nr_spikes = 20,
                      nr_trains = 400,
                      t_stop = None,
                      numpy_seed = None):
    if numpy_seed:
        np.random.seed(numpy_seed)
    temporal_patterns_inter_burst = []
    for x in range(nr_trains):
        # GOTTA DO A .ANY THING HERE!!!
        spike_train = np.array([])
        while not spike_train.any():
            spike_train = np.random.exponential(interval, (nr_spikes)).cumsum()
            if t_stop:
                spike_train = spike_train[spike_train <= t_stop]
        spike_train = np.round(spike_train, decimals = 1)
        temporal_patterns_inter_burst.append(spike_train)

    temporal_patterns_inter_burst = np.array(temporal_patterns_inter_burst)

    return temporal_patterns_inter_burst

def poisson_burst_generator(inter_burst_interval =100,
                            nr_bursts = 20,
                            nr_trains = 400,
                            intra_burst_interval=10,
                            spikes_per_burst = 5,
                            numpy_seed = None,
                            train_t_stop = None,
                            burst_t_stop=None):
    if numpy_seed:
        np.random.seed(numpy_seed)

    # I don't get this fucking .any() thing

    burst_intervals = poisson_generator(interval = inter_burst_interval,
                               nr_spikes = nr_bursts,
                               nr_trains = nr_trains,
                               t_stop = train_t_stop,
                               numpy_seed = numpy_seed)



    new_trains= []
    for train_idx, train in enumerate(burst_intervals):
        bursts = np.array([])
        loop=0
        while not bursts.any():

            bursts = poisson_generator(interval=intra_burst_interval,
                                      nr_spikes=spikes_per_burst,
                                      nr_trains=len(train),
                                      t_stop=burst_t_stop,
                                      numpy_seed=numpy_seed+train_idx+loop)
            loop = loop+1

        bursts = np.round(bursts,decimals=1)
        train = np.repeat(train[:,np.newaxis],np.shape(bursts[1]),axis=1)

        train = np.sort(train)
        new_train = train + bursts
        new_train=np.sort(new_train.flatten())
        new_trains.append(new_train)

    return new_trains

if __name__ == '__main__':
    temporal_patterns_inter_burst = poisson_generator(t_stop=1000, numpy_seed=10000)
    plt.figure()
    plt.eventplot(temporal_patterns_inter_burst)
    plt.title("numpy generated eventplot")

    plt.figure()
    plt.hist(np.hstack(temporal_patterns_inter_burst.flat), bins = int(np.hstack(temporal_patterns_inter_burst.flat).max()))
    plt.title("numpy generated time histogram")
    
    temporal_patterns = poisson_burst_generator(inter_burst_interval=100,
                                               nr_bursts=20,
                                               nr_trains=1000,
                                               intra_burst_interval=100,
                                               spikes_per_burst=3,
                                               numpy_seed=10001,
                                               train_t_stop=1000,
                                               burst_t_stop=None)
    plt.figure()
    plt.eventplot(temporal_patterns)
    plt.title("numpy generated eventplot")

    plt.figure()
    plt.hist(np.hstack(np.array(temporal_patterns).flat), bins = int(np.hstack(np.array(temporal_patterns).flat).max()))
    plt.title("numpy generated time histogram")
    
    isis=[]
    for x in temporal_patterns:
        curr = np.diff(x)
        isis.extend(curr)
    plt.figure()
    isis = np.array(isis)
    plt.hist(isis, bins=int(round(isis.max()/0.1)))
        
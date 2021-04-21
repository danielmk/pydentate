# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 12:49:13 2021

@author: Daniel & barisckuru
"""

import numpy as np
from elephant import spike_train_generation as stg
from neo.core import AnalogSignal
import quantities as pq
from scipy.stats import skewnorm
from skimage.measure import profile_line




def inhom_poiss(modulation_rate=10, max_rate=100, n_cells=400):
    """Generate spike trains from an inhomogeneous poisson process.
    The rate profile is a sine wave with a frequency given my modulation_rate
    and a peak given by max_rate. Both in Hz.
    Returns a ragged array of dtype np.object with lenght n_cells. Contains the
    spike times.
    """
    sampling_interval = 0.0001 * pq.s

    t = np.arange(0, 0.5, sampling_interval.magnitude)

    rate_profile = (np.sin(t*modulation_rate*np.pi*2-np.pi/2) + 1) * max_rate / 2

    rate_profile_as_asig = AnalogSignal(rate_profile,
                                        units=1*pq.Hz,
                                        t_start=0*pq.s,
                                        t_stop=0.5*pq.s,
                                        sampling_period=sampling_interval)

    spike_trains = []
    for x in range(n_cells):
        curr_train = stg.inhomogeneous_poisson_process(rate_profile_as_asig)
        # We have to make sure that there is sufficient space between spikes.
        # If there is not, we move the next spike by 0.1ms
        spike_trains.append(curr_train)

    array_like = np.array([np.around(np.array(x.times)*1000, decimals=1)
                           for x in spike_trains], dtype=np.object)
    for arr_idx in range(array_like.shape[0]):
        bad_idc = np.argwhere(np.diff(array_like[arr_idx]) == 0).flatten()
        bad_idc = bad_idc+1
        while bad_idc.any():
            for bad_idx in bad_idc:
                array_like[arr_idx][bad_idx] = array_like[arr_idx][bad_idx] + 0.1
            bad_idc = np.argwhere(np.diff(array_like[arr_idx]) == 0).flatten()
            bad_idc = bad_idc + 1

    return array_like


def grid_cell_spikes(arena_size=(1, 1), t_span=(0, 500), n_cells=400):
    """Generate spike trains from an inhomogeneous poisson process based on a
    grid cell model with phase precession.
    """
    pass


def _grid_spatial_navigation():
    pass

#Solstad 2006 Grid Model
def _grid_maker(spacing, orientation, pos_peak, arr_size, sizexy, max_rate):
    #define the params from input here, scale the resulting array for maxrate and sperate the xy for size and shift
    arr_size = arr_size #200*200 dp was good enough in terms of resolution
    x, y = pos_peak
    pos_peak = np.array([x,y])
    max_rate = max_rate
    lambda_spacing = spacing*(arr_size/100) #100 required for conversion
    k = (4*np.pi)/(lambda_spacing*np.sqrt(3))
    degrees = orientation
    theta = np.pi*(degrees/180)
    meterx, metery = sizexy
    arrx = meterx*arr_size # *arr_size for defining the 2d array size
    arry = metery*arr_size
    dims = np.array([arrx,arry])
    rate = np.ones(dims)
    #implementation of grid function
    # 3 k values for 3 cos gratings with different angles to generate grid fields
    k1 = ((k/np.sqrt(2))*np.array((np.cos(theta+(np.pi)/12) + np.sin(theta+(np.pi)/12),
          np.cos(theta+(np.pi)/12) - np.sin(theta+(np.pi)/12)))).reshape(2,)
    k2 = ((k/np.sqrt(2))*np.array((np.cos(theta+(5*np.pi)/12) + np.sin(theta+(5*np.pi)/12),
          np.cos(theta+(5*np.pi)/12) - np.sin(theta+(5*np.pi)/12)))).reshape(2,)
    k3 = ((k/np.sqrt(2))*np.array((np.cos(theta+(9*np.pi)/12) + np.sin(theta+(9*np.pi)/12),
          np.cos(theta+(9*np.pi)/12) - np.sin(theta+(9*np.pi)/12)))).reshape(2,)
    
    rate[i,j] = (np.cos(np.dot(k1, curr_dist))+
               np.cos(np.dot(k2, curr_dist))+ np.cos(np.dot(k3, curr_dist)))/3
    rate = max_rate*2/3*(rate+1/2)   # arr is the resulting 2d grid out of 3 gratings
    return rate
    
def _grid_population(n_grid, max_rate, seed, arr_size=200):
    # skewed normal distribution for grid spacings
    np.random.seed(seed)
    median_spc = 43
    spc_max = 100
    skewness = 6  #Negative values are left skewed, positive values are right skewed.
    grid_spc = skewnorm.rvs(a = skewness,loc=spc_max, size=n_grid)  #Skewnorm function
    grid_spc = grid_spc - min(grid_spc)      #Shift the set so the minimum value is equal to zero.
    grid_spc = grid_spc / max(grid_spc)      #Standadize all the vlues between 0 and 1. 
    grid_spc = grid_spc * spc_max         #Multiply the standardized values by the maximum value.
    grid_spc = grid_spc + (median_spc - np.median(grid_spc))
    
    grid_ori = np.random.randint(0, high=60, size=[n_grid,1]) #uniform dist for orientation btw 0-60 degrees
    grid_phase = np.random.randint(0, high=(arr_size-1), size=[n_grid,2]) #uniform dist grid phase
    
    # create a 3d array with grids for n_grid
    rate_grids = np.zeros((arr_size, arr_size, n_grid))#empty array
    for i in range(n_grid):
        x = grid_phase[i][0]
        y = grid_phase[i][1]
        rate = grid_maker(grid_spc[i], grid_ori[i], [x, y], arr_size, [1,1], max_rate)
        rate_grids[:, :, i] = rate
    return rate_grids, grid_spc


def _inhom_poiss(arr, dur_s, poiss_seed=0, dt_s=0.025):
    np.random.seed(poiss_seed)
    n_cells = arr.shape[0]
    spi_arr = np.zeros((n_cells), dtype = np.ndarray)
    for grid_idc in range(n_cells):
        np.random.seed(poiss_seed+grid_idc)
        rate_profile = arr[grid_idc,:]
        asig = AnalogSignal(rate_profile,
                                units=1*pq.Hz,
                                t_start=0*pq.s,
                                t_stop=dur_s*pq.s,
                                sampling_period=dt_s*pq.s,
                                sampling_interval=dt_s*pq.s)
        curr_train = stg.inhomogeneous_poisson_process(asig)
        spi_arr[grid_idc] = np.array(curr_train.times*1000) #time conv to ms
    return spi_arr


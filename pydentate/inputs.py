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
from scipy import stats
import pdb

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d

def sigmoid(x, amplitude, center, spread):
    """
    Compute the sigmoid function for a given input array.

    The sigmoid function is defined as:
    f(x) = amplitude * (1 / (1 + exp((x - center) / spread)))

    Parameters:
    -----------
    x : numpy.ndarray
        The input array for which the sigmoid function will be computed.

    amplitude : float
        The amplitude parameter of the sigmoid function, controlling the vertical scaling.

    center : float
        The center parameter of the sigmoid function, specifying the horizontal shift.

    spread : float
        The spread parameter of the sigmoid function, influencing the slope of the curve.

    Returns:
    --------
    numpy.ndarray
        The result of applying the sigmoid function to the input array `x`.

    Examples:
    ---------
    >>> import numpy as np
    >>> x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    >>> sigmoid(x, amplitude=1.0, center=0.0, spread=1.0)
    array([0.11920292, 0.26894142, 0.5       , 0.73105858, 0.88079708])
    """
    return amplitude * (1 / (1 + np.exp((x - center) / spread)))

def inhomogeneous_poisson_process(
    t_start: float,
    t_stop: float,
    sampling_interval: float,
    rate_profile_frequency: int,
    rate_profile_amplitude: int,
    refractory_period: float
) -> NDArray[np.float64]:
    """
    Generate a spike train from an inhomogeneous Poisson process.

    Parameters
    ----------
    t_start : float
        Start time of the spike train in seconds.
    t_stop : float
        Stop time of the spike train in seconds.
    sampling_interval : float
        Time interval between consecutive time points in seconds.
    rate_profile_frequency : int
        Frequency of the rate profile oscillation in Hz.
    rate_profile_amplitude : int
        Amplitude of the rate profile.
    refractory_period : float | None, optional
        Minimum time between consecutive spikes. If None, no refractory period is enforced.

    Returns
    -------
    spike_times : NDArray[np.float64]
        An array containing spike times generated from the inhomogeneous Poisson process.

    Notes
    -----
    This function generates spike times according to an inhomogeneous Poisson process
    with a specified rate profile. The rate profile is modulated by a sinusoidal function.

    If a refractory period is provided, consecutive spikes will be separated by at least
    the refractory period.

    All units are in seconds.

    Examples
    --------
    >>> spike_train = inhomogeneous_poisson_process(0.0, 10.0, 0.001, 10, 5)
    """

    t = np.arange(t_start, t_stop, sampling_interval)

    rate_profile = rate_profile_amplitude * (np.sin(t * rate_profile_frequency * np.pi * 2 - np.pi / 2) + 1) / 2

    cumulative_rate = np.cumsum(rate_profile) * sampling_interval
    max_cumulative_rate = cumulative_rate[-1]
    n_spikes = np.round(max_cumulative_rate).astype(int)

    random_numbers = np.random.uniform(0, max_cumulative_rate, n_spikes)

    inv_cumulative_rate_func = interp1d(cumulative_rate, t)

    spike_times: NDArray[np.float64] = inv_cumulative_rate_func(random_numbers)
    spike_times.sort()

    if refractory_period is None:
        return spike_times

    thinned_spike_times = np.empty(0, dtype=np.float64)
    previous_spike_time = t_start - refractory_period

    for spike_time in spike_times:
        if spike_time - previous_spike_time > refractory_period:
            thinned_spike_times = np.append(thinned_spike_times, spike_time)
            previous_spike_time = spike_time

    return thinned_spike_times

def homogeneous_poisson_process(
    t_start: float,
    t_stop: float,
    rate: float,
    refractory_period: float
) -> NDArray[np.float64]:
    """
    Generate a spike train from a homogeneous Poisson process.

    Parameters
    ----------
    t_start : float
        Start time of the spike train in seconds.
    t_stop : float
        Stop time of the spike train in seconds.
    rate : float
        Average firing rate of the Poisson process in Hz.
    refractory_period : float | None, optional
        Minimum time between consecutive spikes. If None, no refractory period is enforced.

    Returns
    -------
    spike_times : NDArray[np.float64]
        An array containing spike times generated from the homogeneous Poisson process.

    Notes
    -----
    This function generates spike times from a homogeneous Poisson process with a
    constant firing rate over time. Spike times are randomly drawn from an exponential
    distribution characterized by the provided rate. The spike times are sorted in
    ascending order.

    If a refractory period is provided, consecutive spikes will be separated by at least
    the refractory period.

    All units are in seconds.

    Examples
    --------
    >>> spike_train = homogeneous_poisson_process(0.0, 10.0, 5.0)
    """

    n_spikes = np.random.poisson((t_stop - t_start) * rate)

    spike_times = np.random.uniform(t_start, t_stop, n_spikes)
    spike_times.sort()

    if refractory_period is None:
        return spike_times

    thinned_spike_times = np.empty(0)
    previous_spike_time = t_start - refractory_period

    for spike_time in spike_times:
        if spike_time - previous_spike_time > refractory_period:
            thinned_spike_times = np.append(thinned_spike_times, spike_time)
            previous_spike_time = spike_time

    return thinned_spike_times

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

def gaussian_connectivity_gc_bc(n_pre, n_gc, n_bc, n_syn_gc, n_syn_bc, scale_gc, scale_bc):
    """TODO"""
    pass

def gaussian_connectivity(n_pre, n_post, n_syn=[100,], scale=[1000, 12]):
    """TODO THIS IS A STUB FOR A GENERALIZED VERSION OF gaussian_connectivity_gc_bc
    Choose n_syn postsynaptic cells for each presynaptic cells.
    Clean up this function. It is not Pythonic.
    Possibly vectorize the for loops.

    Parameters
    ----------
    n_pre : int
        Number of presynaptic cells.
    n_post : list
        A list of cell numbers in postsnaptic populations.
    n_syn : int
        Number of synapses from pre to post.
    scale : float
        The scale of the gaussien distribution.

    Returns
    -------
    list : len(n_post)
        Each list entry contains a 2darray of shape (n_pre, n_post).
    """

    center = np.array(n_post) // 2
    gauss = stats.norm(loc=center[0], scale=scale[0])
    pdf = gauss.pdf(np.arange(n_post[0]))
    pdf = pdf/pdf.sum()
    post_idc = np.arange(n_post[0])
    start_idc = np.random.randint(0, n_post[0]-1, size=n_pre)

    out_list = []
    for idx, n_post_pop in enumerate(n_post):
        pre_to_post = []
        # pdb.set_trace()
        for x in start_idc:
            curr_idc = np.concatenate((post_idc[x:n_post_pop], post_idc[0:x]))
            # pdb.set_trace()
            pre_to_post.append(np.random.choice(curr_idc, size=n_syn[idx], replace=False,
                                                p=pdf))
        out_list.append(pre_to_post)
        if idx+1 < len(n_post):
            gauss = stats.norm(loc=center[idx+1], scale=scale[idx+1])
            pdf = gauss.pdf(n_post[idx+1])
            pdf = pdf/pdf.sum()
            post_idc = np.arange(n_post[idx+1])
            start_idc = np.array(((start_idc/n_post_pop)*n_post[idx+1]), dtype=int)


    return np.array(pre_to_post)

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
    
def _grid_population(n_grid, max_rate, seed, arena_size=[1,1], arr_size=200):
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
        rate = _grid_maker(grid_spc[i], grid_ori[i], [x, y], arr_size, arena_size, max_rate)
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


# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 12:49:13 2021

@author: Daniel
"""

import numpy as np
from elephant import spike_train_generation as stg
from neo.core import AnalogSignal
import quantities as pq


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


def _grid_field(spacing):
    pass


def _grid_spatial_navigation():
    pass




# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:46:35 2017

@author: DanielM
"""

from mossycell_cat import MossyCell
from neuron import h, gui
import matplotlib.pyplot as plt
import numpy as np
import quantities as pq
from neo.core import AnalogSignal
import nai.intrinsic_properties as ip
import elephant.statistics as stats

h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\test_granule_cell_Santhakumar2005\\nrnmech.dll")
def _run_current_injections(cell_type, input_current, v_init):
    """DOCSTRING TODO"""
    traces = []
    for x in input_current:
        reload(neuron)
        cell = cell_type()
        soma_v = h.Vector()
        t = h.Vector()
        soma_v.record(cell.soma(0.5)._ref_v)
        stim = h.IClamp(cell.soma(0.5))
        stim.delay = 300
        stim.amp = x
        stim.dur = 500
        t.record(h._ref_t)
        h.cvode.active(0)
        dt = 0.1
        h.steps_per_ms = 1.0/dt
        h.tstop = 1500
        h.finitialize(v_init)
        h.t = -2000
        h.secondorder = 0
        h.dt = 10
        while h.t < -100:
            h.fadvance()
        h.secondorder = 2
        h.t = 0
        h.dt = 0.1
        
        """Setup run control for -100 to 1500"""
        h.frecord_init() # Necessary after changing t to restart the vectors
        while h.t < 1000:
            h.fadvance()
        
        traces.append(soma_v)

    return np.array(traces), t

if __name__ == '__main__':
    v_init = -60
    cell_type = MossyCell
    input_current_passive = np.arange(-0.5,0,0.05)
    input_current_active = np.arange(0,1.1,0.1)

    output_passive, t_vec = _run_current_injections(cell_type, input_current_passive, v_init)
    output_active, t_vec = _run_current_injections(cell_type, input_current_active, v_init)

    # Transform to neo AnalogSignal object
    data_passive = AnalogSignal(np.swapaxes(output_passive,0,1), units = pq.mV, sampling_rate = (1/0.0001) * pq.Hz)
    data_active = AnalogSignal(np.swapaxes(output_active,0,1), units = pq.mV, sampling_rate = (1/0.0001) * pq.Hz)

    # Calculate some statistics
    # Buggy since neo 0.5.2
    """
    minima = ip.minimum_response(data_passive, [300 * pq.ms, 800 * pq.ms])
    steady_state = ip.steady_state(data_passive, [300 * pq.ms, 800 * pq.ms], time = 50 * pq.ms)
    input_resistance_steady = ip.input_R(input_current_passive, steady_state)
    input_resistance_sag = ip.input_R(input_current_passive, minima)
    tau_membrane = ip.tau_membrane(data_passive, [300 * pq.ms, 800 * pq.ms], input_resistance_sag)
    spikes = ip.spike_extraction(data_active)
    avg_freq = []
    for x in spikes:
        avg_freq.append(stats.mean_firing_rate(x, t_start = 300 * pq.ms, t_stop = 800 * pq.ms))
    """
    # Plotting
    for idx in range(np.shape(data_passive)[1]):
        plt.figure()
        curr_sig = data_passive[:,idx]
        plt.plot(curr_sig.times , curr_sig)
        plt.ylim((-100,-55))
        plt.xlabel("time (ms)")
        plt.ylabel("mV")
        fname = "BC_v-init_" + str(v_init) + "passive_" + str(input_current_passive[idx]) + "_uA_.jpg"
        plt.savefig(fname)
        
    for idx in range(np.shape(data_active)[1]):
        plt.figure()
        curr_sig = data_active[:,idx]
        plt.plot(curr_sig.times, curr_sig)
        plt.ylim((-70,55))
        plt.xlabel("time (ms)")
        plt.ylabel("mV")
        fname = "BC_v-init_" + str(v_init) + "active_" + str(input_current_active[idx]) + "_uA_.jpg"
        plt.savefig(fname)
    
        
        
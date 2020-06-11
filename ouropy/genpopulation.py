# -*- coding: utf-8 -*-
"""
This module implements the GenNetwork class, which implements generic network
logic. Also implements Population and GenConnection classes

@author: DanielM
"""

from neuron import h
import random
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import os
import shelve
import scipy.stats as stats


class Population(object):
    """This is the model of a generic population.
    A population is a number of cells of a specific type derived from
    genneuron.GenNeuron. The GenPopulation object keeps track of all
    incoming and outgoing connections. It is recommended to create Populations
    through the GenNetwork.mk_population interface of a network the population
    is part of.

    Attributes
    ----------
    parent_network - gennetwork.GenNetwork or derived instances
        The network the population takes part in
    cell_type - genneuron.GenNeuron class or subclass thereof
        The cell type making up the population
    cells - list of genneuron.GenNeuron instances
        A list of cells that currently exist within the population
    connections - list of Connection objects
        A list of outgoing and incoming connections

    Methods
    -------
    __init__
    make_cells
    get_cell_number
    record_aps
    plot_aps
    write_aps
    current_clamp_rnd
    current_clamp_range
    voltage_recording
    add_connection

    Use cases
    ---------
    >>> nw = GenNetwork()
    >>> nw.mk_population(GranuleCell, 500)
    Create an empty network and create a population of 500 granule cells in the
    network.
    """

    def __init__(self, cell_type=None, n_cells=None, parent_network=None):
        self.parent_network = parent_network
        self.cell_type = cell_type
        self.cells = []
        self.connections = []
        self.VClamps = []
        self.VClamps_i = []
        self.VRecords = []
        if cell_type and n_cells:
            self.make_cells(cell_type, n_cells)
        self.i = 0

    def SEClamp(self, cells, dur1=200, amp1=0, rs=0.001):
        for x in cells:
            clamp = self.cells[x]._SEClamp(dur1=dur1, amp1=amp1, rs=rs)
            self.VClamps.append(clamp)
            curr_vec = h.Vector()
            curr_vec.record(clamp._ref_i)
            self.VClamps_i.append(curr_vec)

    def voltage_recording(self, cells):
        for x in cells:
            record = self.cells[x]._voltage_recording()
            self.VRecords.append(record)

    def make_cells(self, cell_type, n_cells):
        """Create cells of a certain type

        Parameters
        ----------
        cell_type - genneuron.GenNeuron class of subclass thereof
            the type of the cells to be created
        n_cells - numeric
            number of cells to be created

        Returns
        -------
        None

        Use Cases
        ---------
        >>> popul = Population(parent_network = nw)
        >>> popul.make_cells(GranuleCell, 500)
        Create an empty population within nw and then create 500 granule cells
        """

        if hasattr(self, 'cell_type'):
            if self.cell_type != cell_type:
                raise TypeError("cell_type inconsistent with population")
            else:
                self.cell_type = cell_type

        if not hasattr(self, 'cells'):
            self.cells = []

        for _ in range(n_cells):
            self.cells.append(cell_type())

        self.cells = np.array(self.cells, dtype=object)

    def get_cell_number(self):
        """Return the number of cells"""
        return len(self.cells)

    def record_aps(self):
        counters = [cell._AP_counter() for cell in self.cells]
        self.ap_counters = counters
        return counters

    def plot_aps(self, color='k'):
        cells = []
        for x in self.ap_counters:
            # as_numpy() doesn't work on windows 10 ???
            try:
                cells.append(x[0].as_numpy())
            except:
                cells.append(np.array(x[0]))

        # Workaround for matplotlib bug. plt.eventplot throws error when first
        # element empty
        if not np.array(cells[0]).any():
            cells[0] = np.array([0], dtype=float)

        plt.eventplot(cells, linewidth=2, color=color)

    def write_aps(self, directory='', fname=''):
        if not fname:
            time_tup = time.gmtime()
            time_str = time.asctime(time_tup)
            time_str = '_'.join(time_str.split(' '))
            nw_name = self.parent_network.__class__.name
            pop_name = self.cell_type.name
            fname = nw_name + '_' + pop_name + '_' + time_str
            fname = fname.replace(':', '-')
        if not directory:
            directory = os.getcwd()
        if not os.path.isdir(directory):
            os.mkdir(directory)
        path = directory + '\\' + fname + '.npz'
        try:
            ap_list = [x[0].as_numpy() for x in self.ap_counters]
        except:
            ap_list = [np.array(x[0]) for x in self.ap_counters]
        np.savez(path, *ap_list)

    def perc_active_cells(self):
        try:
            # as_numpy doesn't work on windows 10 ???
            timing_arrays = [x[0].as_numpy() for x in self.ap_counters]
        except:
            timing_arrays = [np.array(x[0]) for x in self.ap_counters]
        active_counter = 0
        for x in timing_arrays:
            if x.size != 0:
                active_counter += 1

        return (active_counter / float(self.get_cell_number())) * 100

    def mk_current_clamp(self, cells, amp=0.3, dur=5, delays=3):
        if not hasattr(cells, '__iter__'):
            cells = np.random.choice(self.get_cell_number(), cells,
                                     replace=False)

        if not hasattr(delays, '__iter__'):
            delays = np.array(delays)

        for cell in cells:
            for delay in delays:
                self.cells[cell]._current_clamp_soma(amp=amp, dur=dur,
                                                     delay=delay)

    def current_clamp_rnd(self, n_cells, amp=0.3, dur=5, delay=3):
        """DEPRECATE"""
        chosen_cells = np.random.choice(self.cells, n_cells, replace=False)

        for x in chosen_cells:
            for y in delay:
                x._current_clamp_soma(amp=amp, dur=dur, delay=y)

        return chosen_cells

    def current_clamp_range(self, n_cells, amp=0.3, dur=5, delay=3):
        """DEPRECATE"""
        if type(n_cells) == int:
            n_cells = range(n_cells)

        for cell in n_cells:
            self.cells[cell]._current_clamp_soma(amp=amp, dur=dur, delay=delay)

    """def voltage_recording(self, cell_type):
        rnd_int = random.randint(0, len(self.cells) - 1)
        soma_v_vec = self.cells[rnd_int]._voltage_recording()
        return soma_v_vec"""

    def add_connection(self, conn):
        self.connections.append(conn)

    def get_properties(self):
        """Get the properties of the network"""
        try:
            ap_time_stamps = [x[0].as_numpy() for x in self.ap_counters]
        except:
            ap_time_stamps = [np.array(x[0]) for x in self.ap_counters]
        ap_numbers = [x[1].n for x in self.ap_counters]
        try:
            v_rec = [x.as_numpy() for x in self.VRecords]
            vclamp_i = [x.as_numpy() for x in self.VClamps_i]
        except:
            v_rec = [np.array(x) for x in self.VRecords]
            vclamp_i = [np.array(x) for x in self.VClamps_i]
        properties = {'parent_network': str(self.parent_network),
                      'cell_type': self.cell_type.name,
                      'cell_number': self.get_cell_number(),
                      'connections': [conn.get_properties()
                                      for conn in self.connections],
                      'ap_time_stamps': ap_time_stamps,
                      'ap_number': ap_numbers,
                      'v_records': v_rec,
                      'VClamps_i': vclamp_i}

        properties
        return properties

    def __str__(self):
        return self.cell_type.name + 'Population'

    def __iter__(self):
        return self

    def __getitem__(self, item):
        return self.cells[item]

    def __next__(self):
        if self.i < (len(self.cells)):
            i = self.i
            self.i += 1
            return self.cells[i]
        else:
            self.i = 0
            raise StopIteration()

    def next(self):
        return self.__next__()


# Helpers
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

    return [(x_arr[idx], y_arr[idx]) for idx in range(len(x_arr))]


def euclidian_dist(p1, p2):
    """ p1 and p2 must both be of len 2 where p1 = (x1,y1); p2 = (x2,y2)"""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

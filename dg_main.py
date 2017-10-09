# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 12:36:11 2017

@author: DanielM
"""

"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run
"""

from granulecell import GranuleCell
from mossycell import MossyCell
from basketcell import BasketCell
from hippcell import HippCell
from neuron import h, gui
import matplotlib.pyplot as plt
import ouropy.gennetwork as network
import scalebars as sb
import time

h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\test_granule_cell_Santhakumar2005\\nrnmech.dll")
start_time_offset = time.clock()

#Setup the cells of the network
celltypes = [GranuleCell, MossyCell, BasketCell, HippCell]
cellnums = [10,2,1,1]
my_nw = network.GenNetwork(celltypes, cellnums)

#Setup pp inputs
my_nw.make_connection('pp_stim')
my_nw._make_synapses(h.Exp2Syn, GranuleCell, [4,8], my_nw.connections['pp_stim'])

for x in my_nw.connections['pp_stim'].synapses:
    x.tau1 = 1.5
    x.tau2 = 5.5
    x.e = 0
    x.g = 0.02 #microSiemens

netstim = h.NetStim()
netstim.interval = 100
netstim.number = 1
netstim.start = 10

my_nw.connect_synapses_artificial(netstim, my_nw.connections['pp_stim'], 10, 3,0.1)

#Setup GC to MC connection
my_nw.connect_cells(GranuleCell, MossyCell, [1,5,9,13], divergence = 1,
                        tau1 = 0.5, tau2 = 6.2, e = 0, g_max = 0.0002, thr = 10,
                        delay = 1.5, weight = 0.01, name = "GC->MC")

#Setup GC to BC connection
my_nw.connect_cells(GranuleCell, BasketCell, [1,5,9,13], divergence = 1,
                        tau1 = 0.3, tau2 = 0.6, e = 0, g_max = 0.0047, thr = 10,
                        delay = 0.8, weight = 0.01, name = "GC->BC")

#Setup GC to HC connection
my_nw.connect_cells(GranuleCell, HippCell, [1,4,7,10], divergence = 1,
                    tau1 = 0.3, tau2 = 0.6, e = 0, g_max = 0.5, thr = 10,
                    delay = 1.5, weight = 0.005, name = "GC->HC")

GC_volt, GC_time = my_nw.voltage_recording(GranuleCell)

MC_volt, MC_time = my_nw.voltage_recording(MossyCell)

BC_volt, BC_time = my_nw.voltage_recording(BasketCell)

HC_volt, HC_time = my_nw.voltage_recording(HippCell)

"""Initialization for -2000 to -100"""
h.cvode.active(0)
dt = 0.1
h.steps_per_ms = 1.0/dt
h.tstop = 1500
h.finitialize(-60)
h.t = -2000
h.secondorder = 0
h.dt = 10
while h.t < -100:
    h.fadvance()
    print(h.t)
h.secondorder = 2
h.t = 0
h.dt = 0.1

"""Setup run control for -100 to 1500"""
h.frecord_init() # Necessary after changing t to restart the vectors
while h.t < 150:
    h.fadvance()

exc_time = time.clock()

plt.figure()
plt.plot(GC_time, GC_volt, linewidth = 2.0, color = 'b')
plt.ylim((-80,70))
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)

plt.figure()
plt.plot(MC_time, MC_volt, linewidth = 2.0, color = 'r')
plt.ylim((-80,70))
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)

plt.figure()
plt.plot(BC_time, BC_volt, linewidth = 2.0, color = 'k')
plt.ylim((-80,70))
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)

plt.figure()
plt.plot(HC_time, HC_volt, linewidth = 2.0, color = 'g')
plt.ylim((-80,70))
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)

plt.show()

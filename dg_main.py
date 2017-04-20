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
import network

#h.nrn_load_dll("C:\\nrn\\dentate_gyrus_python_translate\\nrnmech_new.dll")
"""
start_time = time.clock()
n_GCs = 10000

GCs = []
for x in range(n_GCs):
    GCs.append(GranuleCell(h))

end_time = time.clock()
"""
"""
myGC = GranuleCell()

myMC = MossyCell()

myBC = BasketCell()

myHC = HippCell()

GC_volt, GC_time = myGC.somatic_recording()

MC_volt, MC_time = myMC.somatic_recording()

BC_volt, BC_time = myBC.somatic_recording()

HC_volt, HC_time = myHC.somatic_recording()

myGC.simulate()

plt.figure()
plt.plot(GC_time, GC_volt)

plt.figure()
plt.plot(MC_time, MC_volt)

plt.figure()
plt.plot(BC_time, BC_volt)

plt.figure()
plt.plot(HC_time, HC_volt)
"""


celltypes = [GranuleCell]

cellnums = [2]

test_network = network.Network(celltypes, cellnums)

curr_GC = test_network.cells[GranuleCell][0]

curr_GC.pp_stim()

volt, time = curr_GC.somatic_recording()

curr_GC.simulate()

plt.plot(time,volt)

















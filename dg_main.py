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
#import network

#h.nrn_load_dll("C:\\nrn\\dentate_gyrus_python_translate\\nrnmech_new.dll")
"""
myGC = GranuleCell()

myMC = MossyCell()

myBC = BasketCell()

myHC = HippCell()

GC_volt, GC_time = myGC.somatic_recording()

MC_volt, MC_time = myMC.somatic_recording()

BC_volt, BC_time = myBC.somatic_recording()

HC_volt, HC_time = myHC.somatic_recording()

h.tstop = 10000.0
h.run(10000.0)

plt.figure()
plt.plot(GC_time, GC_volt)

plt.figure()
plt.plot(MC_time, MC_volt)

plt.figure()
plt.plot(BC_time, BC_volt)

plt.figure()
plt.plot(HC_time, HC_volt)
"""


#Setup the cells of the network
celltypes = [GranuleCell, MossyCell, BasketCell, HippCell]
cellnums = [10,2,1,1]
my_nw = network.Network(celltypes, cellnums)


#Setup pp inputs
my_nw.make_connection('pp_stim')
my_nw._make_synapses(h.Exp2Sid, GranuleCell, [4,8], my_nw.connections['pp_stim'])
sid = 0

for x in my_nw.connections['pp_stim'].synapses:
    x.tau1 = 1.5
    x.tau2 = 5.5
    x.e = 0
    x.sid = sid
    x.g = 0.02 #microSiemens
    sid += 1

netstim = h.NetStim()
netstim.interval = 100
netstim.number = 1
netstim.start = 500

my_nw.connect_synapses_artificial(netstim, my_nw.connections['pp_stim'], 10, 3, 1)


#Setup GC to MC connection
my_nw.connect_cells(GranuleCell, MossyCell, [1,5,9,13], divergence = 1,
                        tau1 = 0.5, tau2 = 6.2, e = 0, g_max = 0.0002, thr = 10,
                        delay = 1.5, weight = 1, name = "GC->MC")

#Setup GC to BC connection
my_nw.connect_cells(GranuleCell, BasketCell, [1,5,9,13], divergence = 1,
                        tau1 = 0.3, tau2 = 0.6, e = 0, g_max = 0.0047, thr = 10,
                        delay = 0.8, weight = 1, name = "GC->BC")

#Setup GC to HC connection
my_nw.connect_cells(GranuleCell, HippCell, [1,4,7,10], divergence = 1,
                    tau1 = 0.3, tau2 = 0.6, e = 0, g_max = 0.5, thr = 10,
                    delay = 1.5, weight = 1, name = "GC->HC")

GC_volt, GC_time = my_nw.voltage_recording(GranuleCell)

MC_volt, MC_time = my_nw.voltage_recording(MossyCell)

BC_volt, BC_time = my_nw.voltage_recording(BasketCell)

HC_volt, HC_time = my_nw.voltage_recording(HippCell)

#MC_volt, MC_time = my_nw.current_clamp(MossyCell, 0.2, 500, 500)

my_nw.run_network(tstop = 2000)

plt.figure()
plt.plot(GC_time, GC_volt)

plt.plot(MC_time, MC_volt)


plt.show()

"""





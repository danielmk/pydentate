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
celltypes = [GranuleCell, MossyCell]
cellnums = [2, 1]
my_nw = network.Network(celltypes, cellnums)

#Setup pp inputs
my_nw.make_connection('pp_stim')
my_nw.make_synapses(h.Exp2Sid, GranuleCell, [4,8], my_nw.connections['pp_stim'])
sid = 0

for x in my_nw.connections['pp_stim'].synapses:
    x.tau1 = 1.5
    x.tau2 = 5.5
    x.e = 0
    x.sid = sid
    x.g = 0.02 #uS
    sid += 1

netstim = h.NetStim()
netstim.interval = 100
netstim.number = 1
netstim.start = 2000

my_nw.connect_synapses_artificial(netstim, my_nw.connections['pp_stim'], 10, 3, 1)

#Setup GC to MC connection
my_nw.make_connection('GC->MC')
my_nw.make_synapses(h.Exp2Sid, MossyCell, [1,5,9,13], my_nw.connections['GC->MC'])

for x in my_nw.connections['GC->MC'].synapses:
    x.tau1 = 0.5
    x.tau2 = 6.2
    x.e = 0
    x.sid = sid
    x.g = 0.0002
    sid += 1

my_nw.connect_cells(GranuleCell, my_nw.connections['GC->MC'], 1, 10, 1.5 ,1)

GC_volt, GC_time = my_nw.cells[GranuleCell][0].somatic_recording()

#MC_volt, MC_time = my_nw.cells[MossyCell][0].somatic_recording(0,500,500)

myGC = GranuleCell()

new_volt, new_time = myGC.somatic_recording()

h.tstop = 2000

print(h.tstop)

my_nw.run_network()

plt.plot(GC_time, GC_volt)

#plt.plot(MC_time, MC_volt)

"""










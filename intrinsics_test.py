# -*- coding: utf-8 -*-
"""
Created on Mon May 08 12:38:04 2017

@author: DanielM
"""

from granulecell import GranuleCell
from mossycell import MossyCell
from basketcell import BasketCell
from hippcell import HippCell
from neuron import h, gui
import matplotlib.pyplot as plt
import scalebars as sb

myGC = GranuleCell()

myMC = MossyCell()

myBC = BasketCell()

myHC = HippCell()

GC_volt, GC_time = myGC._current_clamp_soma(0.3,500,200)

MC_volt, MC_time = myMC._current_clamp_soma(0.36,500,200)

BC_volt, BC_time = myBC._current_clamp_soma(0.5,500,200)

HC_volt, HC_time = myHC._current_clamp_soma(0.5,500,200)

h.tstop = 800
h.run()

plt.figure()
plt.plot(GC_time, GC_volt)
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)
plt.title("GranuleCell")

plt.figure()
plt.plot(MC_time, MC_volt)
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)
plt.title("MossyCell")

plt.figure()
plt.plot(BC_time, BC_volt)
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)
plt.title("BasketCell")

plt.figure()
plt.plot(HC_time, HC_volt)
plt.axis('off')
sb.add_scalebar(plt.axes(), loc = 1)
plt.title("HippCell")


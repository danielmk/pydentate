# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 10:25:41 2018

@author: daniel
"""
import shelve
import numpy as np
import matplotlib.pyplot as plt

#Home PC
#directory = "C:\\Users\\daniel\\repos\\pyDentate\paradigm_pattern-separation_saves_2018-03-11\\"
#Office PC
#directory = "Y:\\DanielM\\023_Dentate Gyrus Model\\paradigm_spatial-inhibition\\"
#Dropbox
directory = "Z:\\pyDentate\\pyDentateData\\pattern_separation_data_local_input_revised_exp\\seed10006\\scale1000\\net_tunedrev\\"

file_name = "net_tunedrevexp.TunedNetwork_data_paradigm_local-pattern-separation_run_scale_seed_001_1000_10006.pydd"

data = shelve.open(directory + file_name)
# Get to BasketCell Connection
GC_to_BC_targets = data['net_tunedrevexp.TunedNetwork']['populations'][0]['connections'][25]['GranuleCellPopulation to BasketCellPopulation']['pre_cell_targets'].flatten()

# Get to BasketCell Connection
GC_to_HC_targets = data['net_tunedrevexp.TunedNetwork']['populations'][0]['connections'][26]['GranuleCellPopulation to HippCellPopulation']['pre_cell_targets'].flatten()

# Get to BasketCell Connection
BC_to_GC_targets = data['net_tunedrevexp.TunedNetwork']['populations'][0]['connections'][27]['BasketCellPopulation to GranuleCellPopulation']['pre_cell_targets'].flatten()

# Get to HIPPCell Connection
HC_to_GC_targets = data['net_tunedrevexp.TunedNetwork']['populations'][0]['connections'][28]['HippCellPopulation to GranuleCellPopulation']['pre_cell_targets'].flatten()

plt.figure()
plt.hist(GC_to_BC_targets, bins = 24)
plt.xlabel("# BCs")
plt.ylabel("# incoming GC Synapses")
           
plt.figure()
plt.hist(GC_to_HC_targets, bins = 24)
plt.xlabel("# HCs")
plt.ylabel("# incoming GC Synapses")
           
plt.figure()
plt.hist(BC_to_GC_targets, bins = 2000)
plt.xlabel("# GCs")
plt.ylabel("# incoming BC Synapses")

plt.figure()
plt.hist(HC_to_GC_targets,bins = 2000)
plt.xlabel("# GCs")
plt.ylabel("# incoming HC Synapses")
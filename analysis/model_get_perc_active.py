# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 10:25:41 2018

@author: daniel
"""
import shelve
import numpy as np
import matplotlib.pyplot as plt
import os

#Home PC
#directory = "C:\\Users\\daniel\\repos\\pyDentate\paradigm_pattern-separation_saves_2018-03-11\\"
#Office PC
#directory = "Y:\\DanielM\\023_Dentate Gyrus Model\\paradigm_spatial-inhibition\\"
#Dropbox
data_path = "C:\\Users\\DanielM\\Repos\\pyDentate\\pattern_separation_data\\"
file_name = "net_globalrev.TunedNetwork_run_scale_000_500.pydd"
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f and 'norm' in f]
data_files.sort()

data = shelve.open(data_path + file_name)
# Get to BasketCell Connection


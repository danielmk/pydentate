# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h
import numpy as np
import net_globalrev
from burst_generator_inhomogeneous_poisson import inhom_poiss
import os
import argparse
import scipy.stats as stats

# Parse command line inputs
# Command line signature:
# python script -runs n n n -savedir str -scale n -seed n
pr = argparse.ArgumentParser(description='Local pattern separation paradigm')
pr.add_argument('-runs',
                    nargs=3,
                    type=int,
                    help='start stop range for the range of runs',
                    default=[0, 1, 1],
                    dest='runs')
pr.add_argument('-savedir',
                    type=str,
                    help='complete directory where data is saved',
                    default=os.getcwd(),
                    dest='savedir')
pr.add_argument('-scale',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=1000,
                    dest='input_scale')
pr.add_argument('-seed',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=10000,
                    dest='seed')

args = pr.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
input_scale = args.input_scale
seed = args.seed

# Locate a nrnmech.dll file containig the mechanisms required by the network
# adapt path for your own machine.
dll_files = [("C:\\Users\\DanielM\\Repos\\models_dentate\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll"),
              "C:\\Users\\daniel\\Repos\\nrnmech.dll",
             ("C:\\Users\\Holger\\danielm\\models_dentate\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll"),
             ("C:\\Users\\Daniel\\repos\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll")]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + dll_dir)
h.nrn_load_dll(dll_dir)

# Seed the numpy random number generator for replication
np.random.seed(seed)

nw = net_globalrev.TunedNetwork(seed)

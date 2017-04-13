dg_main.py - The main script creating the network
genericcell.py - Defines a generic cell object that encapsualtes the neuron logic of creating sections
connecting the topology, defining geometry and the simulation of a current clamp experiment to test
cell properties
granulecell.py - Defines the specific properties of a granule cell and inherits from genericcell.py
nrnmech.dll - Contains the compiled neuron mechanisms. This file is added to the h object by
default when using from neuron import h, gui
nrnmech_old - A buggy version of the mechanisms, does not give the desired cellular properties

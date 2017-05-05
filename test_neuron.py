# -*- coding: utf-8 -*-
"""
Created on Fri May 05 14:00:13 2017

@author: DanielM
"""

from neuron import h
from neuron import gui

h_attrs = len(dir(h))

print("The h object has " + str(h_attrs) + " attributes.")
print("We have 881 with the proper mechanisms.")
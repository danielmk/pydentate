# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 12:36:11 2017

@author: DanielM
"""

"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run
"""

import numpy

# Set up some simulation parameters
secondorder = 2
tsetp = 0
period = 2
dt = 0.1
tstop = 100     #1500

percentSclerosis = 50
scalingFactor = 500
ngcell = 500
nmcell = 15
nbcell = 5
nhcell = 6
npp = 1
randnet = 0

# Define final network size after sclerosis
nmcell = int(nmcell * ((100 - percentSclerosis) / 100))
nhcell = int(nhcell * ((100 - percentSclerosis ) / 100))

if nmcell == 0:
    nmcell = 1
if nhcell == 0:
    nhcell = 1
    
totalCells = ngcell + nmcell + nbcell + nhcell

#Defining Granule cell
Gcell = [None] * ngcell     #A placeholder list for all granule cells

class GranuleCell(object):
    ndend1 = 4
    ndend2 = 4
    nst = 10
    
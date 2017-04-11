# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 14:45:54 2017

@author: DanielM
"""

import dg_main
#import matplotlib.pyplot as plt

myGC = dg_main.GranuleCell()
volt, time = myGC.somatic_recording()
myGC.simulate()

plt.plot(time, volt)
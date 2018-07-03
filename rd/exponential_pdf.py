# -*- coding: utf-8 -*-
"""
Created on Tue Jul 03 12:13:36 2018

@author: Daniel
"""

import scipy.stats
import matplotlib.pyplot as plt
import numpy as np

pdf_100 = scipy.stats.norm.pdf(np.arange(2000),0,100)
pdf_220 = scipy.stats.expon.pdf(np.arange(2000),0,220)
pdf_400 = scipy.stats.expon.pdf(np.arange(2000),0,400)

cdf_100 = scipy.stats.expon.cdf(np.arange(2000),0,100)
cdf_220 = scipy.stats.expon.cdf(np.arange(2000),0,220)
cdf_400 = scipy.stats.expon.cdf(np.arange(2000),0,400)

plt.figure()
plt.plot(np.arange(2000), pdf_100, color = "r")
plt.plot(np.arange(2000), pdf_220, color = "b")
plt.plot(np.arange(2000), pdf_400, color = "g")

plt.figure()
plt.plot(np.arange(2000), cdf_100, color = "r")
plt.plot(np.arange(2000), cdf_220, color = "b")
plt.plot(np.arange(2000), cdf_400, color = "g")
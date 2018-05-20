# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np

data_path = "C:\\Users\\DanielM\\Repos\\pyDentateData\\frequency_inhibition_data"
onlyfiles = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f))]
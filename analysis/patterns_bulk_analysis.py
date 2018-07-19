# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os

for root, dirs, files in os.walk("R:\\pyDentate\\pyDentateData\\pattern_separation_data_local_input_revised"):
    for name in files:
        if name.endswith('spike_data.npz'):
            data_path = root + '\\'
            # function with input data_path here
            break
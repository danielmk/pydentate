# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import patterns_get_LeutgebMeasure_tresolved_inputs as get_inputs
import patterns_get_LeutgebMeasure_tresolved_outputs as get_outputs
import os

# parent = "R:\\pyDentate\\pyDentateData\\pattern_separation_data_local_input_revised"
parent = "R:\\pyDentate\\pyDentateData\\pattern_separation_data_local_input_revised\\seed10000\\scale100"

for root, dirs, files in os.walk(parent):
    for name in files:
        if name.endswith('spike_data.npz'):
            print(root)
            data_path = root + '\\'
            get_outputs.similarity_measure_leutgeb_output_tresolved_directory(data_path, 1000)
            break
        elif name.startswith('input_patterns') & name.endswith('npz'):
            print(root)
            data_path = root + '\\'
            get_inputs.similarity_measure_leutgeb_inputs_directory(data_path, 1000)
            
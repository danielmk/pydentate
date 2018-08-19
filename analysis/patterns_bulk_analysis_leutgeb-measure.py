# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import patterns_get_LeutgebMeasure_inputs as get_inputs
import patterns_get_LeutgebMeasure_outputs as get_outputs
import os

parent = "Z:\\pyDentate\\pyDentateData\\pattern_separation_data_rate\\net_tunedrev_200\\"
done = 0

for root, dirs, files in os.walk(parent):
    for name in files:
        if os.path.isfile(root + '\\1_leutgeb-measure_matrix_len-bin_6000.txt'):
            done += 1
            break
        elif name.endswith('spike_data.npz'):
            print(root)
            data_path = root + '\\'
            get_outputs.similarity_measure_leutgeb_outputs_directory(data_path, 6000)
            break
        elif name.startswith('input_patterns') & name.endswith('npz'):
            print(root)
            data_path = root + '\\'
            get_inputs.similarity_measure_leutgeb_inputs_directory(data_path, 6000)
            break
print(str(done) + ' files already present')
            
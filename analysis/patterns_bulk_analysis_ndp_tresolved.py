# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import patterns_get_NDP_inputs_tresolved as get_inputs
import patterns_get_NDP_outputs_tresolved as get_outputs
import os

parent = "Z:\\pyDentate\\pyDentateData\\pattern_separation_data_local_30Hz_input\\"
done = 0

for root, dirs, files in os.walk(parent):
    for name in files:
        if os.path.isfile(root + '\\1_NDP_matrix_333.txt'):
            print(root + " that one is DONE")
            done += 1
            break
        elif name.endswith('spike_data.npz'):
            print(root)
            data_path = root + '\\'
            get_outputs.similarity_measure_NDP_directory(data_path, 1000)
            break
        elif name.startswith('input_patterns') & name.endswith('npz'):
            print(root)
            data_path = root + '\\'
            get_inputs.similarity_measure_NDP_directory(data_path, 1000)
            break
print(str(done) + ' files already present')

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import model_get_spike_stats_active
import os

parent = "Z:\\pyDentate\\pyDentateData\\pattern_separation_data_local_input_revised\\"
done = 0

for root, dirs, files in os.walk(parent):
    for name in files:
        if os.path.isfile(root + '\\1_perc_active_cells.txt'):
            print(root + " that one is DONE")
            done += 1
            break
        elif name.endswith('spike_data.npz'):
            print(root)
            data_path = root + '\\'
            model_get_spike_stats_active.get_spike_stats(data_path)
            break
"""
        elif name.startswith('input_patterns') & name.endswith('npz'):
            print(root)
            data_path = root + '\\'
            get_inputs.similarity_measure_NDP_directory(data_path)
            break
"""
print(str(done) + ' files already present')

# -*- coding: utf-8 -*-
"""
@author: DanielM
"""

from ouropy import gennetwork
import numpy as np
from granulecell import GranuleCell
from mossycell_cat import MossyCell
from basketcell import BasketCell
from hippcell import HippCell


class TunedNetwork(gennetwork.GenNetwork):
    """This network class is identical to net_tunedrev.TunedNetwork except that
    the synaptic weight from all interneurons are set to zero."""
    name = "TunedNetwork"

    def __init__(self, seed=None, temporal_patterns=np.array([]),
                 spatial_patterns_gcs=np.array([]),
                 spatial_patterns_bcs=np.array([])):
        self.init_params = locals()
        self.init_params['self'] = str(self.init_params['self'])
        # Setup cells
        self.mk_population(GranuleCell, 2000)
        self.mk_population(MossyCell, 60)
        self.mk_population(BasketCell, 24)
        self.mk_population(HippCell, 24)

        # Set seed for reproducibility
        if seed:
            self.set_numpy_seed(seed)

        # Setup recordings
        self.populations[0].record_aps()
        self.populations[1].record_aps()
        self.populations[2].record_aps()
        self.populations[3].record_aps()

        temporal_patterns = np.array(temporal_patterns)
        if (type(spatial_patterns_gcs) == np.ndarray and
           type(temporal_patterns) == np.ndarray):

            for pa in range(len(spatial_patterns_gcs)):
                # PP -> GC
                gennetwork.PerforantPathPoissonTmgsyn(self.populations[0],
                                                      temporal_patterns[pa],
                                                      spatial_patterns_gcs[pa],
                                                      'midd', 10, 0, 1, 0, 0,
                                                      1*10**(-3))

        if (type(spatial_patterns_bcs) == np.ndarray and
           type(temporal_patterns) == np.ndarray):

            for pa in range(len(spatial_patterns_bcs)):
                # PP -> BC
                gennetwork.PerforantPathPoissonTmgsyn(self.populations[2],
                                                      temporal_patterns[pa],
                                                      spatial_patterns_bcs[pa],
                                                      'ddend', 6.3, 0, 1, 0, 0,
                                                      1*10**(-3))

        # GC -> MC
        gennetwork.tmgsynConnection(self.populations[0], self.populations[1],
                                    12, 'proxd', 1, 7.6, 500, 0.1, 0, 0, 10,
                                    1.5, 0.2*10**(-2) * 10)

        # GC -> BC
        # Weight x4, target_pool = 2
        gennetwork.tmgsynConnection(self.populations[0], self.populations[2],
                                    8, 'proxd', 1, 8.7, 500, 0.1, 0, 0, 10,
                                    0.8, 2.5*10**(-2))

        # GC -> HC
        # Divergence x4; Weight doubled; Connected randomly.
        gennetwork.tmgsynConnection(self.populations[0], self.populations[3],
                                    24, 'proxd', 1, 8.7, 500, 0.1, 0, 0, 10,
                                    1.5, 2.5*10**(-2))

        # MC -> MC
        gennetwork.tmgsynConnection(self.populations[1], self. populations[1],
                                    24, 'proxd', 3, 2.2, 0, 1, 0, 0, 10,
                                    2, 0.5*10**(-3))

        # MC -> BC
        gennetwork.tmgsynConnection(self.populations[1], self.populations[2],
                                    12, 'proxd', 1, 2, 0, 1, 0, 0, 10,
                                    3, 0.3*10**(-3))

        # MC -> HC
        gennetwork.tmgsynConnection(self.populations[1], self.populations[3],
                                    20, 'midd', 2, 6.2, 0, 1, 0, 0, 10,
                                    3, 0.2*10**(-3))

        # BC -> GC
        # Nr. synapses x3; Weight *1/4;tau from 5.5 to 20 (Hefft & Jonas, 2005)
        gennetwork.tmgsynConnection(self.populations[2], self.populations[0],
                                    560, 'soma', 400, 20, 0, 1, 0, -70, 10,
                                    0.85, 0)

        # We reseed here to make sure that following connections are consistent
        # between this and net_globalrev which has a global target pool for
        # BC->GC.
        if seed:
            self.set_numpy_seed(seed+1)

        # BC -> MC
        gennetwork.tmgsynConnection(self.populations[2], self.populations[1],
                                    28, 'proxd', 3, 3.3, 0, 1, 0, -70, 10,
                                    1.5, 0)

        # BC -> BC
        gennetwork.tmgsynConnection(self.populations[2], self.populations[2],
                                    12, 'proxd', 2, 1.8, 0, 1, 0, -70, 10,
                                    0.8, 0)

        # HC -> GC
        # Weight x10; Nr synapses x4; tau from 6 to 20 (Hefft & Jonas, 2005)
        gennetwork.tmgsynConnection(self.populations[3], self.populations[0],
                                    2000, 'dd', 640, 20, 0, 1, 0, -70, 10,
                                    3.8, 0)

        # HC -> MC
        gennetwork.tmgsynConnection(self.populations[3], self.populations[1],
                                    60, ['mid1d', 'mid2d'], 4, 6, 0, 1, 0, -70,
                                    10, 1, 0)

        # HC -> BC
        gennetwork.tmgsynConnection(self.populations[3], self.populations[2],
                                    24, 'ddend', 4, 5.8, 0, 1, 0, -70, 10,
                                    1.6, 0)

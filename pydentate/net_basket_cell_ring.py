# -*- coding: utf-8 -*-
"""
This module implements the class TunedNetwork.
TunedNetwork creates a ring network as defined in Santhakumar et al. 2005
with some changes as in Yim et al. 2015.
See StandardNetwork docstring for details.
Created on Tue Nov 28 13:01:38 2017

@author: DanielM
"""

from ouropy import gennetwork
import numpy as np
from . import granulecell
from . import mossycell_cat
from . import  basketcell
from . import hippcell
from .inputs import sigmoid

GranuleCell = granulecell.GranuleCell
MossyCell = mossycell_cat.MossyCell
BasketCell = basketcell.BasketCell
HippCell = hippcell.HippCell

class BasketCellRing(gennetwork.GenNetwork):
    """A ring network of basket cells. Represents about 300um of dorsal rat
    hippocampus. 
    Huskoo et al. 2015: 4000 PVIR cells in entire unilateral hippocampus.
    By rough approximation about 1cm from medial to temporal tip of hippocampus.
    Cells per 300um: (4000 / 1cm) * 0.03cm = 120 PVIR Cells
    
    """
    name = "BasketCellRing"

    def __init__(self, seed, temporal_patterns, spatial_patterns_bcs, rec_matrix, gap_matrix,
                 rec_weight=7.6e-3, n_bcs=120, gap_resistance=6e2, gap_delay=0.0, gap_junctions=True, pp_bc_weight=1e-3, pv_reversal=-70):
        self.init_params = locals()
        self.init_params['self'] = str(self.init_params['self'])

        # Setup cells
        self.mk_population(BasketCell, n_bcs)

        # Set seed for reproducibility
        if seed:
            self.set_numpy_seed(seed)

        # Setup recordings
        self.populations[0].record_aps()

        temporal_patterns = np.array(temporal_patterns, dtype=object)

        # weights
        # feedforward inhibition

        for pa in range(len(spatial_patterns_bcs)):
            # PP -> BC
            gennetwork.PerforantPathPoissonTmgsyn(self.populations[0],
                                                  temporal_patterns[pa],
                                                  spatial_patterns_bcs[pa],
                                                  'ddend', 6.3, 0, 1, 0, 0,
                                                  pp_bc_weight)



        # BC -> BC
        gennetwork.tmgsynConnectionMatrix(self.populations[0], self.populations[0],
                                    rec_matrix, 'proxd', 1.8, 0, 1, 0, pv_reversal, 10,
                                    0.8, rec_weight)
        
        if gap_junctions:
            gennetwork.GapJunctionConnectionMatrix(self.populations[0], self.populations[0], gap_matrix, 'proxd', gap_resistance, gap_delay)

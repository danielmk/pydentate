# -*- coding: utf-8 -*-
"""
This script provides a rough template for a script that can fit the synpatic
conductances measured at a simulated granule cell to the experimentally
measured conductances.

It needs to perform several steps:
1. Load the experimental data
2. Setup a granule cell and a synaptic process with a parameterset w
3. Simulate the synaptic conductance at the granule cell
4. Calculate a distance metric between stimulated and measured conductance
5. Optimize w to minimize the distance metric

There are many ways to do this. In the script I make SUGGESTIONS for the
components the script should contain. I also cannot guarantee that I thought
of everything.

@author: Daniel
"""

import scipy.io
import scipy.optimize
from granulecell import GranuleCell

"""
LOAD THE EXPERIMENTAL DATA

The data from Martin is in .mat files.
scipy.io.loadmat provides an interface to load .mat files. Once loaded the 
excitatory conductance should be saved under a key called 'ge'.
"""


"""
SETUP/SIMULATE THE GRANULE CELL AND SYNAPTIC PROCESS

Setup and simulation can happen within the same function. The parameterset that
defines the outcome of the simulation is (tau_1, tau_2, tau_facil, tau_rec).
Those parameters define the synaptic properties. The intrinsic properties are
set by GranuleCell.
"""

def simulate(tau_1, tau_2, tau_facil, tau_rec):
    # Create GranuleCell
    # Attach tmgexp2syn to GranuleCell
    # Setup recording from GranuleCell soma with h.SEClamp()
    # Run Simulation
    # Return output arrays
    pass

"""
CALCULATE A DISTANCE METRIC

The distance metric (I call is loss function, which comes from deep learning)
is supposed to express the difference between the experimental ge and the
simulated ge. There are infinitely many ways to calculate the "distance" 
between two arrays. A commonly used one loss function is the mean squared 
error, with which it is probably worth starting.

IMPORTANT: Remember that the loss function should be invariant to the absolute
amounts of conductance. That means that both y1 and y2 will have to be
normalized in some way to account for absolute effects.
ALSO: To compare two arrays meaningfully, they must have the same number of
datapoints, so sampling rate matters!
"""

def loss(y1,tau_1, tau_2, tau_facil, tau_rec, method="mse"):
    pass

"""
OPTIMIZATION

scipy.optimize should get the job done here. I am not sure yet if the
computation time will be a problem, since we essentially have to setup and
destroy a neuron model at every iteration of the optimization. But this has to
be seen in practice.
IMPORTANT: Initial values will matter for the optimization.
"""

    
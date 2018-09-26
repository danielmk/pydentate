# ouropy

ouropy is a generic wrapper of [NEURON](https://www.neuron.yale.edu/neuron/), intended to simplify network modeling. When working in native NEURON, the logical unit is the section.
A section is a cable with physical dimensions and biophysical mechanisms. Building a network means to create sections, connect them to form neurons, populate them with mechanism and
finally connect them synaptically to finish the network. While NEURON makes this rather easy, we still need to keep track of all our sections. The goal of ouropy is to move from the
logic of sections to the logic of neurons, populations and finally entire networks.

# Dependencies

This project mainly depends on a running version of [NEURON](https://www.neuron.yale.edu/neuron/).
Python needs to be able to import the neuron module. Therefor it also needs to have the binaries of NEURON in
its search path. For details refer to [NEURONs official documentation](https://www.neuron.yale.edu/neuron/docs).

# License

This project is published under a [GPL v2](http://www.gnu.org/licenses/gpl-2.0.html)
For a full copy of the license refer to the file GPLv2-LICENSE.md

# Author
Daniel Müller, Institute of Experimental Epileptology and Cognition Research

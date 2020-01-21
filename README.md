# pyDentate

pyDentate is a biophysically realistic computational model of the dentate gyrus, a hippocampal brain region associated with memory formation and a computation called pattern separation. The properties of the this model are based on a [dentate model] (http://www.opensourcebrain.org/projects/dentate) from the Soltesz Laboratory. We made changes based on new literature and our own experimental findings. Furthermore, we introduced enhancements to study pattern separation.

# Getting started
Follow these steps to get started with pyDentate
<ol>
<li>Install <a href="https://www.anaconda.com/distribution">Anaconda</a></li>
<li>Install <a href="https://www.neuron.yale.edu/neuron">NEURON</a>
  <p>There are many ways to install NEURON. I prefer the <a href="https://anaconda.org/conda-forge/neuron">conda-forge</a> distribution<blockquote>
        <p>conda install -c conda-forge/label/cf201901 neuron</p>
    </blockquote></p>
</li>
<li>Install elephant
  <p><blockquote>pip install elephant</blockquote></li></p>
<li><a href="https://www.neuron.yale.edu/neuron/download/compile_mswin">Compile the NEURON mechanisms</a> in /mechs_7-6</li>
<li>Download the pyDentateeLife2020 repository and unpack</li>
<li>Open paradigm_pattern_separation.py and add the path to your compiled mechanisms to dll_files variable</li>
<li>Run a paradigm_ file</li>
</ol>

# Usage
To run an existing model you simple have to execute a paradigm_ file after following the setup steps above. paradigm_ files are scripts that take care of everything from creating networks, simulating the networks and saving the output. Networks in turn

# References
pyDentate builds on a computational model from this paper: Santhakumar, V., Aradi, I., & Soltesz, I. (2005). Role of mossy fiber sprouting and mossy cell loss in hyperexcitability: a network model of the dentate gyrus incorporating cell types and axonal topography. Journal of Neurophysiology, 93(1), 437–453. https://doi.org/10.1152/jn.00777.2004

# Author

Daniel Müller-Komorowska - [Institute of Experimental Epileptology and Cognition Research](https://eecr-bonn.de/)

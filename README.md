# pydentate

pydentate is a biophysically realistic computational model of the dentate gyrus, a hippocampal brain region associated with memory formation and a computation called pattern separation.  We made changes based on new literature and our own experimental findings. Furthermore, we introduced enhancements to study pattern separation.

# Installation
<ol>
<li><p>Open a terminal and clone the pydentate repository</p><blockquote>git clone https://github.com/danielmk/pydentate.git</blockquote></li>
<li><p>Inside the cloned repository pip install</p><blockquote>cd pydentate<br />pip install -e .</blockquote>ON WINDOWS: You will need to <a href="https://www.neuron.yale.edu/neuron/download">install neuron manually.</a></li>

  <li><p>The installation is now complete and you can get started by running the example script</p><blockquote>python paradigm_pattern_separation_baseline.py</blockquote>This will take a long time. If no errors are raised your pydentate is working.</li>
</ol>
If you encounter problems with running pydentate or have questions feel free to contact me (muellerkomorowska@protonmail.com
 or https://twitter.com/scidanm).

# References
pyDentate builds on a computational model from this paper: Santhakumar, V., Aradi, I., & Soltesz, I. (2005). Role of mossy fiber sprouting and mossy cell loss in hyperexcitability: a network model of the dentate gyrus incorporating cell types and axonal topography. Journal of Neurophysiology, 93(1), 437–453. https://doi.org/10.1152/jn.00777.2004
Their model can be found [here](http://www.opensourcebrain.org/projects/dentate).

We decribe pydentate in our eLife paper: Braganza, O., Müller-Komorowska, D., Kelly, T., & Beck, H. (2020). Quantitative properties of a feedback circuit predict frequency-dependent pattern separation. ELife, 813188. https://doi.org/10.7554/eLife.53148
You can find a separated repository that is dedicated to reproducing the paper results [here](https://github.com/danielmk/pyDentateeLife2020).

# Authors

Daniel Müller-Komorowska - [Institute of Experimental Epileptology and Cognition Research](https://eecr-bonn.de/)

Barış Can Kuru (synaptic_fitting add on) - [Institute of Experimental Epileptology and Cognition Research](https://eecr-bonn.de/)

# pyDentate

pyDentate is a biophysically realistic computational model of the dentate gyrus, a hippocampal brain region associated with memory formation and a computation called pattern separation. The properties of the this model are based on a [dentate model] (http://www.opensourcebrain.org/projects/dentate) from the Soltesz Laboratory. We made changes based on new literature and our own experimental findings. Furthermore, we introduced enhancements to study pattern separation.

# Usage
To get started with pyDentate you will need to resolve some depndencies

To use pyDentate for yourself or reproduce our findings, you will need to install [NEURON](https://www.neuron.yale.edu/neuron/) and your python environment needs to be aware
of the \nrn\lib\python\neuron package that is installed with NEURON and the NEURON .dlls need to be available to python. Finally, you will need to compile the NEURON meachnisms
that come with pyDentate in pyDenate/mechs on your machine with NEURONs mknrndll tool (nrnivmodl on linux) which will give you a nrnmech.dll file. The paradigm scripts run networks and in these
scripts you need to give the path to that nrnmech.dll file.

All published results were produced with NEURON7.4 and Python 2.7
The Model also runs in Python 3 and Neuron 7.5 and 7.6 but for unknown reasons the exact spike times do not replicate.
Refer to pyDentate\output_examples to see outputs generated with different versions.

If you experience problems with running pydentate contact danielmuellermsc@gmail.com.

# Author

Daniel Müller-Komorowska - [Institute of Experimental Epileptology and Cognition Research](https://eecr-bonn.de/)

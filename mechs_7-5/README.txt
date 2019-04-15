These mechanisms are included with pyDentate for purposes of replicating the model.
Most of these mechanisms were taken from the Santhakumar (2005) model.
To run pyDentate you have to compile the mechanism on the machine you want to run the
model on. To do so you have to install NEURON and use its mknrndll program on the
folder containing the mechanisms. If succesfull, this will create a nrnmech.dll.
Most parts of pyDentate need to be made aware where that nrnmech.dll is located on
your machine.
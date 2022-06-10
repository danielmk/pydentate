import os

from . import net_tunedrev, inputs, input_generator

dirname = os.path.dirname(__file__)

linux_precompiled = os.path.join(dirname, 'x86_64', 'libnrnmech.so')
windows_precompiled = os.path.join(dirname, 'win64', 'nrnmech.dll')

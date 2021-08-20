import os

from . import net_tunedrev, inputs, input_generator

dirname = os.path.dirname(__file__)

linux_precompiled = os.path.join(os.path.split(dirname)[0], 'x86_64', 'libnrnmech.so')
windows_precompiled = os.path.join(os.path.split(dirname)[0], 'win64', 'nrnmech.dll')

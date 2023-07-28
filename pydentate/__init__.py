import os

dirname = os.path.dirname(__file__)

linux_precompiled = os.path.join(dirname, "x86_64", ".libs", "libnrnmech.so")
windows_precompiled = os.path.join(dirname, "win64", "./libs", "nrnmech.dll")

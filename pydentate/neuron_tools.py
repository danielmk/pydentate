from neuron import h, gui
from pydentate import windows_precompiled, linux_precompiled
import platform

def load_compiled_mechanisms(path='precompiled'):
    """Loads precompiled mechanisms in pyDentate unless
    path defines the full path to a compiled mechanism file."""
    if path != 'precompiled':
        h.nrn_load_dll(path)
    else:
        if platform.system() == 'Windows':
            h.nrn_load_dll(windows_precompiled)
        else:
            h.nrn_load_dll(linux_precompiled)

def run_neuron_simulator(warmup=2000, dt_warmup=10, dt_sim=0.1, t_start=0,
                        t_stop=600, v_init=-60):
    h.cvode.active(0)
    dt = 0.1
    h.steps_per_ms = 1.0/dt
    h.finitialize(v_init)
    h.t = -warmup
    h.secondorder = 0
    h.dt = dt_warmup
    while h.t < -100:
        h.fadvance()

    h.secondorder = 2
    h.t = 0
    h.dt = dt_sim

    """Setup run control for -100 to 1500"""
    h.frecord_init()  # Necessary after changing t to restart the vectors
    while h.t < t_stop:
        h.fadvance()

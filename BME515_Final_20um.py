# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 12:38:02 2022

Created by Eric Musselman January 29, 2021
@author: Winnie Lu, Philjae Chang, Celina Zhou

"""

# %% import NEURON library
from neuron import h
import matplotlib.pyplot as plt

# other imports
import numpy as np
import csv

# %% load run controls
h.load_file('stdrun.hoc')

# %% MODEL SPECIFICATION
# time params =========================================================================
h.tstop = 10  # [ms]: simulation time
h.dt = 0.001  # [ms]: timestep

# cell params =========================================================================
h.celsius = 37
# [mV]: Vm @ rest for initializing membrane potential at start of simulation
Vo = -80
n_nodes = 81  # []: (int) number of sections, make this an odd number
D = 20  # [um]: fiber diameter
inl = 100*D  # [um]: internodal length
rhoa = 54.7  # [Ohm]: axoplasmic/axial resistivity
cm = 2.5  # [uF/cm**2]
L = 1.5  # [um]
nseg = 1  # []: (int)
g = 1/2000  # [S/cm**2]

# stim params =========================================================================
delay = 0.5               # [ms]: start time of stim
dur = 1                # [ms]: pulse width of (monopolar) stim
# [mA]: amplitude of (intracellular) stim object -- but we are applying this extracellular (negative cathodic, positive anodic)


vol = []
# change this filename for different diameters
filename = "20um_fiber.txt"
f = open(filename, "r")
lines = f.readlines()
vol = []
for x in lines:
    if x[0:1] != "%":
        vol.append(float(x.split(" ")[-1].strip("\n")))
f.close()


# %% MODEL INITIALIZATION
# define nodes for cell =================================================================
nodes = [h.Section(name=f'node[{i}]') for i in range(n_nodes)]

# insert extracellular/mechanisms part of the circuit ===================================
# connect the nodes =====================================================================
for node_ind, node in enumerate(nodes):
    node.nseg = nseg
    node.diam = 0.6*D
    node.L = L
    node.Ra = rhoa*((L+inl)/L)
    node.cm = cm
    node.insert('sweeney')
    node.insert('extracellular')
    for seg in node:
        seg.extracellular.e = 0
    if node_ind > 0:
        node.connect(nodes[node_ind-1](1))

# %% INSTRUMENTATION - STIMULATION/RECORDING

vol_mem = [h.Vector().record(sec(0.5)._ref_v) for sec in nodes]
tvec = h.Vector().record(h._ref_t)

# %% SIMULATION CONTROL
# compute extracellular potentials from point current source (call this from my_advance to update at each timestep)


def update_field():
    phi_e = []
    for node_ind, node in enumerate(nodes):
        phi_e.append(vol[node_ind])
        node(0.5).e_extracellular = 10000 * phi_e[node_ind] * \
            35 * np.sin(2 * np.pi * 5000 * h.t * 0.001)

# time integrate with constant time step - this just defines method, called by proc advance() below


def my_advance():
    update_field()
    h.fadvance()


h.finitialize(Vo)

# this is somewhat of a "hack" to change the default run procedure in HOC ==================
h(r"""
proc advance() {
    nrnpython("my_advance()")
}""")

# run until tstop ==========================================================================
h.continuerun(h.tstop)

# %% DATA POST PROCESSING / OUTPUT
# plot things ==============================================================================
print(vol_mem)

plt.figure(num=1, clear=False)
plt.plot(tvec, vol_mem[40])  # plot membrane potential at node 10
plt.xlabel('time (ms)')
plt.ylabel('Vm (mV)')
plt.title("gammaCore Activation at 20 um")
plt.show()

print('============ DONE ============')

print('============ DONE ============')

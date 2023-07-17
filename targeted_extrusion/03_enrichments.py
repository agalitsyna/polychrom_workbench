import pickle
import os
import time
import numpy as np
import polychrom

from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file

import simtk.openmm
import os
import shutil


import warnings
import h5py
import glob

from loadutils import *
from simutils import *
import matplotlib.pyplot as plt

from sys import argv

input = argv[1]
output = argv[2]

if not output.endswith("pdf"):
    raise Exception(f"Output file format unknown: {output}")

# -------defining parameters----------
#  -- basic loop extrusion parameters
myfile = h5py.File(input, mode="r")

N = myfile.attrs["N"]
LEFNum = myfile.attrs["LEFNum"]
LEFpositions = myfile["positions"]

mtx = np.zeros([5000, 5000])

for locs in LEFpositions[:, :, :]:
    for pos1, pos2 in locs:
        if np.isfinite(pos1) and np.isfinite(pos2):
            mtx[int(pos1 % 5000), int(pos2 % 5000)] += 1
            mtx[int(pos2 % 5000), int(pos1 % 5000)] += 1

plt.figure(figsize=(13, 13))
plt.imshow(mtx, cmap="coolwarm", vmin=-0, vmax=0.5)
plt.colorbar()
plt.savefig(f"{output}")

# Aggregated and pileup plots:
from mirnylib import numutils

positions = [1000, 2000, 3000, 4000]
aggcoef = 10

for window in [400]:
    pile = np.zeros([window, window])
    pile_agg = np.zeros([window // aggcoef, window // aggcoef])

    for i in range(len(positions)):
        position = positions[i]
        fragment = mtx[
            position - window // 2 : position + window // 2,
            position - window // 2 : position + window // 2,
        ]
        pile += fragment

        pile_agg += numutils.coarsegrain(fragment, aggcoef)

    plt.figure(figsize=(13, 13))
    plt.imshow(pile, cmap="coolwarm", vmin=-0, vmax=np.nanpercentile(pile, 97))
    plt.xticks(
        [0, 10, 20, 30, 40],
        [
            f"{x:.0f} Kb"
            for x in [
                -window // aggcoef,
                -0.5 * window // aggcoef,
                0,
                0.5 * window // aggcoef,
                window // aggcoef,
            ]
        ],
    )
    plt.colorbar()
    plt.savefig(f"{output}.pile.{window}.pdf")

    plt.figure(figsize=(13, 13))
    plt.imshow(pile_agg, cmap="coolwarm", vmin=-0, vmax=np.nanpercentile(pile_agg, 97))
    plt.xticks(
        [0, 10, 20, 30, 40],
        [
            f"{x:.0f} Kb"
            for x in [
                -window // aggcoef,
                -0.5 * window // aggcoef,
                0,
                0.5 * window // aggcoef,
                window // aggcoef,
            ]
        ],
    )
    plt.colorbar()
    plt.savefig(f"{output}.pile_agg.{window}.pdf")


#### Cohesin profile:
coh_profile = np.zeros(5000)
for locs in LEFpositions[:, :, :]:
    for pos1, pos2 in locs:
        if np.isfinite(pos1) and np.isfinite(pos2):
            coh_profile[int(pos1 % 5000)] += 1
            coh_profile[int(pos2 % 5000)] += 1

positions = [1000, 2000, 3000, 4000]
window = 400  # full window size
aggcoef = 10

profile_pile = np.zeros(window)
profile_pile_agg = np.zeros(window // aggcoef)

for i in range(len(positions)):
    position = positions[i]
    fragment = coh_profile[position - window // 2 : position + window // 2]
    profile_pile += fragment

    profile_pile_agg += numutils.coarsegrain(fragment, aggcoef)

plt.figure(figsize=(10, 5))
plt.plot(profile_pile_agg / np.nanmean(profile_pile_agg), color="k")
plt.plot(
    np.arange(0, 40, 1 / aggcoef),
    profile_pile / np.nanmean(profile_pile),
    color="k",
    alpha=0.2,
)
plt.xticks(
    [0, 10, 20, 30, 40],
    [
        f"{x:.0f} Kb"
        for x in [
            -window // aggcoef,
            -0.5 * window // aggcoef,
            0,
            0.5 * window // aggcoef,
            window // aggcoef,
        ]
    ],
)
plt.savefig(f"{output}.cohesin-profile.{window}.pdf")

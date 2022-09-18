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
import time

from loadutils import *
from simutils import *

### Reading input and output names:
from sys import argv
infile = argv[1]
folder = argv[2]  # for output

### Load input 1D trajectory:
myfile = h5py.File(infile, mode="r")

N = myfile.attrs["N"]
LEFNum = myfile.attrs["LEFNum"]
LEFpositions = myfile["positions"]
Nframes = LEFpositions.shape[0]

steps = 200  # MD steps per step of cohesin, same as in Ed's paper
stiff = 1
dens = 0.1
box = (N / dens) ** 0.33  # density = 0.1.
data = grow_cubic(N, int(box) - 2)  # creates a compact conformation
block = 0  # starting block

# new parameters because some things changed
saveEveryBlocks = 10  # save every 10 blocks (saving every block is now too much almost)
restartSimulationEveryBlocks = 100  # in the Ed's paper this was set to 200

# parameters for smc bonds
smcBondWiggleDist = 0.2
smcBondDist = 0.5

# assertions for easy managing code below
assert (Nframes % restartSimulationEveryBlocks) == 0
assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0

savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
simInitsTotal = (Nframes) // restartSimulationEveryBlocks

### Starting simulations and measuring time for simulations
start_time = time.time()

milker = bondUpdater(LEFpositions)

reporter = HDF5Reporter(
    folder=folder,
    max_data_length=100,
    overwrite=True,
    blocks_only=False,
    check_exists=False,
)

for iteration in range(simInitsTotal):

    # simulation parameters are defined below
    a = Simulation(
        platform="cuda",
        integrator="variableLangevin",
        error_tol=0.01,
        GPU="2",
        collision_rate=0.03,
        N=len(data),
        reporters=[reporter],
        PBCbox=[box, box, box],
        precision="mixed",
    )  # timestep not necessary for variableLangevin

    ############################## New code ##############################
    a.set_data(data)  # loads a polymer, puts a center of mass at zero

    a.add_force(
        forcekits.polymer_chains(
            a,
            chains=[(0, None, 0)],
            # By default the library assumes you have one polymer chain
            # If you want to make it a ring, or more than one chain, use self.setChains
            # self.setChains([(0,50,1),(50,None,0)]) will set a 50-monomer ring and a chain from monomer 50 to the end
            bond_force_func=forces.harmonic_bonds,
            bond_force_kwargs={
                "bondLength": 1.0,
                "bondWiggleDistance": 0.1,  # Bond distance will fluctuate +- 0.05 on average
            },
            angle_force_func=forces.angle_force,
            angle_force_kwargs={
                "k": 1.5
                # K is more or less arbitrary, k=4 corresponds to presistence length of 4,
                # k=1.5 is recommended to make polymer realistically flexible; k=8 is very stiff
            },
            nonbonded_force_func=forces.polynomial_repulsive,
            nonbonded_force_kwargs={
                "trunc": 1.5,  # this will let chains cross sometimes
                "radiusMult": 1.05,  # this is from old code
                #'trunc':10.0, # this will resolve chain crossings and will not let chain cross anymore
            },
            except_bonds=True,
        )
    )

    # ------------ initializing milker; adding bonds ---------
    # copied from addBond
    kbond = a.kbondScalingFactor / (smcBondWiggleDist**2)
    bondDist = smcBondDist * a.length_scale

    activeParams = {"length": bondDist, "k": kbond}
    inactiveParams = {"length": bondDist, "k": 0}
    milker.setParams(activeParams, inactiveParams)

    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(
        bondForce=a.force_dict["harmonic_bonds"], blocks=restartSimulationEveryBlocks
    )

    # If your simulation does not start, consider using energy minimization below
    if iteration == 0:
        a.local_energy_minimization()
    else:
        a._apply_forces()

    for i in range(restartSimulationEveryBlocks):
        if i % saveEveryBlocks == (saveEveryBlocks - 1):
            a.do_block(steps=steps)
        else:
            a.integrator.step(
                steps
            )  # do steps without getting the positions from the GPU (faster)
        if i < restartSimulationEveryBlocks - 1:
            curBonds, pastBonds = milker.step(
                a.context
            )  # this updates bonds. You can do something with bonds here
    data = a.get_data()  # save data and step, and delete the simulation
    del a

    reporter.blocks_only = True  # Write output hdf5-files only for blocks

    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

reporter.dump_data()

end_time = time.time()
total_time = end_time - start_time


print(f"{total_time} sec, {total_time/60:.2f} min")

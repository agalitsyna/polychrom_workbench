import pickle
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

import click

### Reading input and output names:
@click.command()
@click.option("-v", "--verbose", count=True)
@click.option(
    "-f",
    "--force/--no-force",
    default=False,
    help="Forced rewriting output folder. ",
)
@click.option(
    "--infile",
    "-i",
    required=True,
    help=" Input 1D trajectory in hdf5 format. ",
    type=click.File('r'),
)
@click.option(
    "--output-folder",
    "-o",
    required=True,
    help=" Output folder for writing 3D simualtions output. ",
    type=str,
)
@click.option(
    "--gpu",
    "-G",
    default=0,
    show_default=True,
    help="GPU index",
    type=int)
@click.option(
    "--steps",
    "-s",
    default=200,
    show_default=True,
    help="MD steps per step of cohesin, 200 is recommended default based on Banigan et al. 2022",
    type=int,
)
@click.option(
    "--dens",
    "-d",
    default=0.1,
    show_default=True,
    type=np.float32,
    help="Initial density of the polymer " "(will be constructed by grow_cubic)",
)
@click.option(
    "--error-tol",
    "-e",
    default=0.01,
    show_default=True,
    type=np.float32,
    help="Error tolerance in the simulations ",
)
@click.option(
    "--collision-rate",
    "-c",
    default=0.03,
    show_default=True,
    type=np.float32,
    help="Collision rate in simulations ",
)
@click.option(
    "--smc-bond-wiggle-dist",
    "-w",
    default=0.2,
    show_default=True,
    type=np.float32,
)
@click.option(
    "--smc-bond-dist",
    "-b",
    default=0.5,
    show_default=True,
    type=np.float32,
)
def run_3D_simulations(
    infile,
    output_folder,
    gpu,
    force,
    verbose,
    steps,
    dens,
    error_tol,
    smc_bond_wiggle_dist,
    smc_bond_dist,
    collision_rate,
):
    ### Load input 1D trajectory:
    myfile = h5py.File(infile.name, mode="r")

    N = myfile.attrs["N"]
    LEFpositions = myfile["positions"]
    Nframes = LEFpositions.shape[0]

    box = (N / dens) ** 0.33  # density = 0.1.
    data = grow_cubic(N, int(box) - 2)  # creates a compact conformation

    # new parameters because some things changed
    saveEveryBlocks = (
        10  # save every 10 blocks (saving every block is now too much almost)
    )
    restartSimulationEveryBlocks = 100  # in the Ed's paper this was set to 200

    # parameters for smc bonds
    smcBondWiggleDist = smc_bond_wiggle_dist  # 0.2
    smcBondDist = smc_bond_dist  # 0.5

    # assertions for easy managing code below
    assert (Nframes % restartSimulationEveryBlocks) == 0
    assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0

    simInitsTotal = (Nframes) // restartSimulationEveryBlocks

    ### Starting simulations and measuring time for simulations
    start_time = time.time()

    milker = bondUpdater(LEFpositions)

    if force:
        if os.path.exists(output_folder):
            if verbose:
                print(f" {output_folder} exists, rewriting in --force mode.")
            os.remove(output_folder)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    reporter = HDF5Reporter(
        folder=output_folder,
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
            error_tol=error_tol,
            GPU=f"{gpu}",
            collision_rate=collision_rate,
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
            bondForce=a.force_dict["harmonic_bonds"],
            blocks=restartSimulationEveryBlocks,
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

    if verbose:
        print(f"{total_time} sec, {total_time/60:.2f} min")


if __name__ == "__main__":
    run_3D_simulations()

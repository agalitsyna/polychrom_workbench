import sys
sys.path.append("../")
from loadutils import *
from simutils import *

import pyximport; pyximport.install(setup_args={"script_args":["--compiler=unix"], "include_dirs":np.get_include()}, reload_support=True)
from smcTranslocator_MovingBarrier import smcTranslocatorDirectional

import logging
logging.basicConfig(encoding="utf-8", level=logging.INFO)

import click

import numpy as np
import random
import joblib
import h5py
import tqdm
import os

@click.command()
@click.option(
    '-v',
    '--verbose',
    count=True)
@click.option(
    "--lifetime",
    "-L",
    default=100,
    show_default=True,
    help="lifetime",
    type=int)
@click.option(
    "--separation",
    "-S",
    default=200,
    help="separation of non-targeted LEFs",
    type=int)
@click.option(
    "--system-size",
    "-N",
    default=200,
    show_default=True,
    help="size of the single system ",
    type=int
)
@click.option(
    "--n-repeats",
    "-R",
    default=200,
    show_default=True,
    help="how many times to repeat the system",
    type=int,
)
@click.option(
    "--trajectory-length",
    "-D",
    required=True,
    help="trajectory length (duration of simulations in number of steps)",
    type=int,
)
@click.option(
    "--step",
    default=50,
    show_default=True,
    help="step length of simulations (how frequently to slice the sims)",
    type=int,
)
@click.option(
    "--outfile",
    "-o",
    "--output",
    required=True,
    help="Output file",
    type=click.File('wb'),
)
@click.option(
    "--pol-kin",
    default=0.01,
    show_default=True,
    help="PolII initiation rate, Default range: 0.00025-0.002",
    type=np.float32,
)
@click.option(
    "--pol-kter",
    default=1.0,
    show_default=True,
    help="Termination rate of PolII, Default range: 0.002-1.0. "
         "Normally transcription is initiation limited, so make sure that pol_kin < pol_kter",
    type=np.float32,
)
@click.option(
    "--pol-step-prob",
    default=0.1,
    show_default=True,
    help="The speed of PolII as a fraction of the speed of cohesin. ",
         #"Step probability of PolII, usually speed set to 80bp/s, "
         #"cohesin speed=1500bp/s, according to Banigan et al. 2022",
    type=np.float32,
)
@click.option(
    "--pol-loading",
    type=str,
    default="35,96,115,176",
    help="List of transcription start sites in a single system (comma-separated list of integers)"
)
@click.option(
    "--pol-kin-loading",
    type=str,
    default="0.5,0.25,1.5,1",
    help="List of kin at transcription start sites in a single system (comma-separated list of integers)."
         "Values are the fractions out of full pol_kin for each gene in pol_loading"
)
@click.option(
    "--pol-termination",
    type=str,
    default="65,66,145,146",
    help="List of transcription termination sites in a single system (comma-separated list of integers) "
         "Make sure that genes don't overlap!"
)
@click.option(
    "--pol-dissociation",
    default=0,
    show_default=True,
    help="Polymerase dissociation probability. ",
    type=np.float32,
)
@click.option(
    "--lef-unstall-rate",
    default=1/(2*100), #0.005,
    show_default=True,
    help="Unstalling rate of LEF. Usually is about once per two typical LEF lifetimes",
    type=np.float32,
)
@click.option(
    "--lef-permeability",
    default=0,
    show_default=True,
    help="controls amount of LEF-LEF bypassing we have",
    type=np.float32,
)
@click.option(
    "--collision-factor",
    default=1.0,
    show_default=True,
    help="factor by which head-on collision with RNAP changes lifetime of cohesin",
    type=np.float32,
)
@click.option(
    "--codir-collision-factor",
    default=1.0,
    show_default=True,
    help="factor for co-directional collisions",
    type=np.float32,
)
@click.option(
    "--pol-permeability-left",
    default=0.,
    show_default=True,
    help="Set permeability of PolII to cohesin. L=0 means PolII is impermeable.",
    type=np.float32,
)
@click.option(
    "--pol-permeability-right",
    default=0.,
    show_default=True,
    help="Set permeability of PolII to cohesin coming from the right",
    type=np.float32,
)


def run_moving_barriers(
    lifetime,
    separation,
    system_size,
    n_repeats,
    trajectory_length,
    step,
    outfile,
    verbose,
    # Simulations-specific parameters:
    pol_kin,
    pol_kter,
    pol_step_prob,
    pol_loading,
    pol_kin_loading,
    pol_termination,
    pol_dissociation,
    lef_unstall_rate,
    lef_permeability,
    collision_factor,
    codir_collision_factor,
    pol_permeability_left,
    pol_permeability_right,
):

    # Verifying parameters:
    if verbose==0:
        tqdm_disable = True
    else:
        tqdm_disable = False

    ### Define full system size:
    N = system_size * n_repeats

    logging.info(f"Loading parameters... System size: {N}")

    pol_loading = np.array(pol_loading.split(',')).astype(int)
    pol_kin_loading = np.array(pol_kin_loading.split(',')).astype(np.float32)
    pol_termination = np.array(pol_termination.split(',')).astype(int)
    assert pol_loading.shape==pol_kin_loading.shape, "Length of Pol loading sites should be the same as annotation of their kinteics"
    assert pol_loading.shape==pol_termination.shape, "Length of Pol loading sites should be the same as length of termination sites"

    pol_loading_full = pol_loading.copy()
    pol_termination_full = pol_termination.copy()

    # Extend pol loading/termination/Kin to full size of the system (with repeats):
    for i in range(n_repeats-1):
        pol_loading_full = np.concatenate([pol_loading_full, pol_loading + i*system_size])
        pol_termination_full = np.concatenate([pol_termination_full, pol_termination + i*system_size])
    pol_kin_loading_full = np.tile(pol_kin_loading, n_repeats)

    ### Define arrays defining behavior of the system:
    # Describe polymerase loading and unloading kinetics:
    pol_kinArray = np.zeros(N, dtype=np.double)
    pol_kterArray = np.zeros(N, dtype=np.double)

    for i, pos in enumerate(pol_loading_full):
        pol_kinArray[pos] = pol_kin * pol_kin_loading_full[i]

    PolSpeedArray = pol_step_prob * np.ones(N, dtype=np.double)

    # Describe polymerase stalling
    stalProbPol = np.zeros(N, dtype=np.double)
    # Note from Ed/Aafke: The rate of PolII stalling. Once stalled, PolII unloads with a rate 'pol_kter'.
    # One can choose a single stall site, or a range of stall sites. If this array is set to zero everywhere,
    # PolII will stall at the defined termination sites.
    unstallArray = np.zeros(N, dtype=np.double)  # array for unstalling in of Pol II in gene

    # Modify stalling, unstalling at TSS/TES.
    # Here we assume that polymerase stops at TES.
    # The following two parameters are hard-coded for Dicty:
    stall_in_gene = 0.
    unstall_in_gene = 1.
    for pos_start, pos_end in zip(pol_loading_full, pol_termination_full):
        # gene on the forward strand:
        if pos_end > pos_start:
            stalProbPol[pos_start:pos_end] = stall_in_gene
            unstallArray[pos_start:pos_end] = unstall_in_gene
            pol_kterArray[pos_end] = pol_kter
        # gene on the reverse strand:
        else:
            stalProbPol[pos_end:pos_start+1] = stall_in_gene
            unstallArray[pos_end:pos_start+1] = unstall_in_gene
            pol_kterArray[pos_end] = pol_kter


    # Cohesin loading properties:
    birthArray = np.zeros(N, dtype=np.double) + 0.1  # cohesin loading, relative units
    deathArray = np.zeros(N, dtype=np.double) + 1. / (lifetime)

    # Stalling of cohesin (e.g., CTCF, which we don't add for Dicty):
    stallLeftArray = np.zeros(N, dtype=np.double)
    stallRightArray = np.zeros(N, dtype=np.double)
    stallDeathArray = np.zeros(N, dtype=np.double) + 1. / (lifetime) # set unbinding rate during stall all the same

    unstallLEFArray = np.ones(N, dtype=np.double) * lef_unstall_rate
    # Usually set a single unstall rate. Non-stall sites can have unstall>0, since there won't be stalls hereanyway.

    pauseArray = np.ones(N, dtype=np.double)
    # The speed of LEFS at each position. Default is an array of ones, which means that LEFS go on maximum speed everywhere.

    shrinkPauseArray = np.zeros(N, dtype=np.double)  # speed at which loops shrink

    pauseArrayPol = np.zeros(N, dtype=np.double) + PolSpeedArray  # speed of PolII at each lattice site

    # Describe collisions of cohesin and RNA Pol:
    collisionLife = lifetime * collision_factor # change upon head-on collision
    coDirCollisionLife = lifetime * codir_collision_factor # change upon co-directional collision
    collisionDeathArray = np.zeros(N, dtype=np.double) + 1. / collisionLife
    coDirCollisionDeathArray = np.zeros(N, dtype=np.double) + 1. / coDirCollisionLife

    # Describe permeability of polymerase to cohesin:
    permLeftArray = np.zeros(N, dtype=np.double)
    permRightArray = np.zeros(N, dtype=np.double)
    permLeftArray += pol_permeability_left
    permRightArray += pol_permeability_right

    LEFNum = N // separation

    logging.info("Starting 1D simulations...")
    SMCTran = smcTranslocatorDirectional(
                      birthArray,
                      deathArray,
                      stallLeftArray,
                      stallRightArray,
                      unstallLEFArray,
                      pauseArray,
                      stallDeathArray,
                      LEFNum,
                      pol_kinArray,
                      pol_kterArray,
                      pauseArrayPol,
                      shrinkPauseArray,
                      pol_loading_full,
                      pol_termination_full,
                      stalProbPol,
                      unstallArray,
                      PolPermL=permLeftArray,
                      PolPermR=permRightArray,
                      collisionFalloffProb=collisionDeathArray,
                      coDirCollisionFalloffProb=coDirCollisionDeathArray,
                      poldissoc=pol_dissociation,
                      LefPerm=lef_permeability,
                      strongCTCF=0 #strongCTCFstall
                      )

    ### Saving trajectories
    with h5py.File(outfile.name, mode="w") as myfile:
        try:
            dset = myfile.create_dataset(
                "positions",
                shape=(trajectory_length, LEFNum, 2),
                dtype=np.float32,
                compression="gzip",
            )
            dset_pol = myfile.create_dataset("PolOccupancy",
                                             shape=(trajectory_length, N),
                                             dtype=np.int32,
                                             compression="gzip")

            for i in tqdm.tqdm(range(trajectory_length), disable=tqdm_disable ):
                cur = []
                cur_pol = []

                for _ in range(step):
                    SMCTran.steps(1)

                LEFs = SMCTran.getSMCs()
                cur.append(np.array(LEFs).T)

                Pols = SMCTran.getPolOccupied()
                cur_pol.append(np.array(Pols))

                cur = np.array(cur)
                dset[i : i + 1] = cur

                cur_pol = np.array(cur_pol)
                dset_pol[i : i + 1] = cur_pol

            myfile.attrs["N"] = N
            myfile.attrs["n_lef"] = LEFNum

            myfile.attrs.update({
                "lifetime": lifetime,
                "separation": separation,
                "system_size": system_size,
                "n_repeats": n_repeats,
                "trajectory_length": trajectory_length,
                "outfile": outfile.name,
                "n_lef": LEFNum,
                "pol_kin": pol_kin,
                "pol_kter": pol_kter,
                "pol_step_prob": pol_step_prob,
                "pol_loading": pol_loading,
                "pol_kin_loading": pol_kin_loading,
                "pol_termination": pol_termination,
                "pol_dissociation": pol_dissociation,
                "lef_unstall_rate": lef_unstall_rate,
                "lef_permeability": lef_permeability,
                "collision_factor": collision_factor,
                "codir_collision_factor": codir_collision_factor,
                "pol_permeability_left": pol_permeability_left,
                "pol_permeability_right": pol_permeability_right,
            })

        except Exception as e:
            if os.path.exists(outfile.name):
                os.remove(outfile.name)

            raise Exception(e)


if __name__ == "__main__":
    run_moving_barriers()

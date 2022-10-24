import os.path

from polychrom.hdf5_format import list_URIs, load_URI

from polychrom.contactmaps import (
    monomerResolutionContactMap,
    monomerResolutionContactMapSubchains,
    binnedContactMap,
)
from multiprocessing import Pool
import glob
from polychrom import polymerutils
import numpy as np
import time

import logging
logging.basicConfig(encoding="utf-8", level=logging.INFO)

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
    "--input-folder",
    "-i",
    required=True,
    help=" Input 3D trajectory folder. ",
    type=str,
)
@click.option(
    "--output",
    "-o",
    required=True,
    help=" Output file for writing 3D simualtions output. Will be in .npy format",
    type=str,
)
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
    "--expected-trajectory-length",
    "-D",
    default=None,
    show_default=True,
    help="Expected length of the trajectory (for assertion)",
    type=int,
)

def get_2D_contactmap(
    input_folder,
    output,
    force,
    verbose,
    system_size,
    n_repeats,
    expected_trajectory_length
):

    ### Starting simulations and measuring time for simulations
    start_time = time.time()

    output_file = output if output.endswith("npy") else output+".npy"

    if os.path.exists(output_file):
        if not force:
            raise Exception(f"{output_file} exists, remove it or use --force option.")
        elif verbose:
            logging.info(f"{output_file} exists, rewriting with --force option.")


    N = system_size * n_repeats
    starts = list(range(0, N, system_size))
    URIs = list_URIs(input_folder)

    logging.info(f"Trajectory length is ... {len(URIs)}")

    if not expected_trajectory_length is None:
        assert len(URIs)==expected_trajectory_length, "The simulations are not completed yet... re-submit"

    hmap = monomerResolutionContactMapSubchains(
        filenames=URIs,
        mapStarts=starts,
        mapN=system_size,
        cutoff=5,  # 2 corresponds to 50-100 nm, 1 monomer
        n=8,
        loadFunction=lambda x: load_URI(x)["pos"],
    )

    np.save(output_file, hmap)


    end_time = time.time()
    total_time = end_time - start_time

    if verbose:
        logging.info(f"Finished in: {total_time} sec, {total_time/60:.2f} min")


if __name__ == "__main__":
    get_2D_contactmap()

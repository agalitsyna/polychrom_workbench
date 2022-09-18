from loadutils import *
from simutils import *

import logging

logging.basicConfig(encoding="utf-8", level=logging.INFO)
logging.info("Starting 1D simulations...")

import tqdm

import click


@click.command()
@click.option(
    '-v',
    '--verbose',
    count=True)
@click.option(
    "--lifetime",
    "-L",
    default=100,
    help="lifetime",
    type=int)
@click.option(
    "--separation",
    "-S",
    default=None, # recommended default is 200
    help="separation of non-targeted LEFs",
    type=int)
@click.option(
    "--separation-total",
    default=None,
    help="separation of all LEFs, targeted plus non-targeted."
         " Cannot be provided together with regular --separation",
    type=int)
@click.option(
    "--n-lef-targeted",
    "-T",
    default=None,
    help="number of targeted LEFs", type=int
)
@click.option(
    "--enrichment",
    "-E",
    default=None,
    help="relative enrichmetn at loading sites."
         " Keep in mind that n_lef/n_lef_t = f.s./enrichment where f.s. is "
         "fountain separation (mean distance between loading sites). "
         " Cannot be provided with number of targeted LEFs --n-lef-targeted",
    type=np.float32,
)
@click.option(
    "--loading-sites",
    type=str,
    default="",
    help="List of loading sites in each system (comma-separated list of integers)",
)
@click.option(
    "--tads",
    type=str,
    default="",
    help="List of tads in each system (comma-separated list of integers)"
)
@click.option(
    "--ctcf-capture-probability",
    default=0.25,
    help="CTCF capture probability",
    type=np.float32,
)
@click.option(
    "--ctcf-release-probability",
    default=0.005,
    help="CTCF release probability",
    type=np.float32,
)
@click.option(
    "--system-size",
    "-N",
    default=6000,
    help="size of the single system ",
    type=int
)
@click.option(
    "--n-repeats",
    "-R",
    default=20,
    help="how many times to repeat the system",
    type=int,
)
@click.option(
    "--trajectory-length",
    "-D",
    default=None,
    help="trajectory length (duration of simulations in number of steps)",
    type=int,
)
@click.option(
    "--step",
    default=10,
    help="step length of simulations (how frequently to slice the sims)",
    type=int,
)
@click.option(
    "--outfile",
    "-o",
    "--output",
    required=True,
    help="Output file")


def run_targeted_loading_extrusion(
    lifetime,
    separation,
    separation_total,
    n_lef_targeted,
    enrichment,
    loading_sites,
    tads,
    ctcf_capture_probability,
    ctcf_release_probability,
    system_size,
    n_repeats,
    trajectory_length,
    step,
    outfile,
    verbose
):

    # Verifying parameters:
    if verbose==0:
        tqdm_disable = True
    else:
        tqdm_disable = False

    if enrichment is None and n_lef_targeted is None:
            raise ValueError("Either enrichment or n_lef_targeted should be provided")
    if not enrichment is None and not n_lef_targeted is None:
        raise ValueError("Enrichment and number of targeted LEFs cannot be provided together")

    if separation is None and separation_total is None:
        raise ValueError("Either separation or total separation of LEFs should be provided")
    if not separation is None and not separation_total is None:
        raise ValueError("Both separation and total separation cannot be provided together")

    if len(loading_sites)>0:
        loading_sites = [int(x) for x in loading_sites.split(',')]
    else:
        loading_sites = []

    if len(tads)>0:
        tads = [int(x) for x in tads.split(',')]
    else:
        tads = []

    # Calculating characteristics of the simulations:
    L = system_size * n_repeats  # System size
    L_t = len(loading_sites) * n_repeats  # Number of beads with enriched loading
    fountain_separation = system_size / len(loading_sites)

    if not separation is None: # Setting separation for only non-targeted LEFs
        n_lef = L // separation  # Number of regular LEFs
        if n_lef==0:
            n_lef = 1
        if not enrichment is None:
            n_lef_targeted = int( enrichment * L_t * n_lef // L )
        if not n_lef_targeted is None:
            enrichment = n_lef_targeted * L / (L_t * n_lef)

        n_lef_total = n_lef_targeted + n_lef
        if n_lef_total==0:
            raise ValueError("Total number of LEFs cannot be zero")
        separation_total = L / n_lef_total

    elif not separation_total is None:
        n_lef_total = int( L // separation_total )
        if n_lef_total == 0:
            raise ValueError("Total number of LEFs cannot be zero")
        if not n_lef_targeted is None:
            n_lef = n_lef_total - n_lef_targeted
            if n_lef == 0:
                n_lef = 1
            enrichment = fountain_separation * n_lef_targeted / n_lef
        elif not enrichment is None:
            n_lef = int( n_lef_total // ( 1 + enrichment/fountain_separation ) )
            if n_lef == 0:
                n_lef = 1
            n_lef_targeted = n_lef_total - n_lef
        separation = L / n_lef # Define separation for non-targeted LEFs

    if n_lef_targeted==0:
        separation_t = 100500
        separation_t_v1 = 100500
    else:
        separation_t = L_t / n_lef_targeted  # Derived separation of the targeted LEFs
        separation_t_v1 = L / n_lef_targeted

    logging.info(
        f"""
Run targeted loading simulations with parameters:
    lifetime, {lifetime}
    separation, {separation}
    separation_total, {separation_total}
    n_lef_targeted, {n_lef_targeted}
    enrichment, {enrichment}
    loading_sites, {loading_sites}
    tads, {tads}
    ctcf_capture_probability, {ctcf_capture_probability}
    ctcf_release_probability, {ctcf_release_probability}
    system_size, {system_size}
    n_repeats, {n_repeats}
    trajectory_length, {trajectory_length}
    step, {step}
Writing to outfile: {outfile}
Derived parameters:
    separation_t, {separation_t}
    separation_t_v1, {separation_t_v1}
    fountain_separation, {fountain_separation}
    n_lef, {n_lef}
    n_lef_total, {n_lef_total}
"""
    )

    # Setting parameters:
    ctcfLeftRelease = {}
    ctcfRightRelease = {}
    ctcfLeftCapture = {}
    ctcfRightCapture = {}

    for i in range(n_repeats):
        for tad in tads:
            pos = i * system_size + tad
            ctcfLeftCapture[pos] = ctcf_capture_probability
            ctcfLeftRelease[pos] = ctcf_release_probability
            ctcfRightCapture[pos] = ctcf_capture_probability
            ctcfRightRelease[pos] = ctcf_release_probability

    ### Regular loading:
    args = {}
    args["N"] = L
    args["LIFETIME"] = lifetime
    args["LIFETIME_STALLED"] = lifetime  # no change in lifetime when stalled
    args["ctcfRelease"] = {-1: ctcfLeftRelease, 1: ctcfRightRelease}
    args["ctcfCapture"] = {-1: ctcfLeftCapture, 1: ctcfRightCapture}

    probsLoading = np.ones(L) / L
    for i in range(n_repeats):
        for site in loading_sites:
            pos = i * system_size + site
            probsLoading[pos] = 0  # no loading of regular LEFs on targeted sites
    probsLoading = probsLoading / np.sum(probsLoading)
    args["probsLoading"] = probsLoading

    # Tergeted_loading
    args_targeted = {}
    args_targeted["N"] = L
    args_targeted["LIFETIME"] = lifetime
    args_targeted["LIFETIME_STALLED"] = lifetime  # no change in lifetime when stalled
    args_targeted["ctcfRelease"] = {-1: ctcfLeftRelease, 1: ctcfRightRelease}
    args_targeted["ctcfCapture"] = {-1: ctcfLeftCapture, 1: ctcfRightCapture}

    probsLoading_targeted = np.zeros(L)
    for i in range(n_repeats):
        for site in loading_sites:
            pos = i * system_size + site
            probsLoading_targeted[pos] = 1
    probsLoading_targeted = probsLoading_targeted / np.sum(probsLoading_targeted)
    args_targeted["probsLoading"] = probsLoading_targeted

    occupied = np.zeros(L)
    # Formally setting ends of the system:
    occupied[0] = 1
    occupied[-1] = 1
    occupied[1] = 1
    occupied[-2] = 1

    # Load targeted cohesins:
    cohesins_targeted = []
    for i in range(n_lef_targeted):
        loadOneWithProbs(cohesins_targeted, occupied, args_targeted)

    # Load all the rest:
    cohesins = []
    for i in range(n_lef):
        loadOneWithProbs(cohesins, occupied, args)

    with h5py.File(outfile, mode="w") as myfile:
        dset = myfile.create_dataset(
            "positions",
            shape=(trajectory_length, n_lef + n_lef_targeted, 2),
            dtype=np.float32,
            compression="gzip",
        )

        for i in tqdm.tqdm(range(trajectory_length), disable=tqdm_disable ):
            cur = []

            for _ in range(step):
                translocate(cohesins, occupied, args)
                translocate(cohesins_targeted, occupied, args_targeted)

            positions = [
                (cohesin.left.pos, cohesin.right.pos)
                for cohesin in cohesins + cohesins_targeted
            ]

            cur.append(positions)
            cur = np.array(cur)
            dset[i : i + 1] = cur

        myfile.attrs["N"] = L
        myfile.attrs["n_lef"] = n_lef

        myfile.attrs.update({
            "lifetime": lifetime,
            "separation": separation,
            "separation_total": separation_total,
            "n_lef_targeted": n_lef_targeted,
            "enrichment": enrichment,
            "loading_sites": ','.join([str(x) for x in loading_sites]),
            "tads": ','.join([str(x) for x in tads]),
            "ctcf_release_probability": ctcf_release_probability,
            "ctcf_capture_probability": ctcf_capture_probability,
            "system_size": system_size,
            "n_repeats": n_repeats,
            "trajectory_length": trajectory_length,
            "outfile": outfile,
            "separation_t": separation_t,
            "separation_t_v1": separation_t_v1,
            "fountain_separation": fountain_separation,
            "n_lef": n_lef,
            "n_lef_total": n_lef_total,
        })


if __name__ == "__main__":
    run_targeted_loading_extrusion()

# loading_sites = [1000, 2000, 3000, 4000, 5000]

# for separation in [1000, 500, 250, 200, 100, 50,]:
#    for lifetime in [100, 50, 200, 300, 500, 1000]:
#        for n_lef_targeted in [2, 5, 10, 50]:

# run_targeted_loading_extrusion(
#    N1=6000,
#    M=20,
#    LIFETIME=lifetime,
#    SEPARATION=separation,
#    n_lef_targeted=n_lef_targeted,
#    trajectoryLength=101000,
#    loadingSites=loading_sites,
#    tads=[x for y in loading_sites for x in [y-150, y+150]],
#    stepLength=10,
#    outfile=f"trajectory/LEFPositions.targeted.{lifetime}.{separation}.{n_lef_targeted}.{CTCF_capture_prob}.h5"
# )

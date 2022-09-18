from polychrom.hdf5_format import list_URIs, load_URI

from polychrom.contactmaps import monomerResolutionContactMap, monomerResolutionContactMapSubchains, binnedContactMap
from multiprocessing import Pool 
import glob
from polychrom import polymerutils
import numpy as np


from mirnylib import numutils


from sys import argv
input_folder = argv[1]
output_file = argv[2] +(".npy" if not argv[2].endswith("npy") else "")


N1 = 2500 # size of one system
M = 20
N = N1 * M 

starts = list(range(0, N, N1))

URIs = list_URIs(input_folder)
print(len(URIs))

hmap = monomerResolutionContactMapSubchains(
    filenames=URIs,
    mapStarts=starts,
    mapN=N1,
    cutoff=5,# 2 corresponds to 50-100 nm, 1 monomer
    n=8,
    loadFunction=lambda x:load_URI(x)["pos"])

np.save(output_file, hmap)

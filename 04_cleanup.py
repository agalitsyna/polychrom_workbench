import h5py
import glob
import os
from sys import argv

if len(argv)>1:
    dry_run = True
else:
    dry_run = False

files_traj1d = glob.glob("data/traj1d/*.h5py")

for infile in files_traj1d:
    try:
        myfile = h5py.File(infile, mode="r")

        N = myfile.attrs["N"]

        myfile.close()

        print("Ok... ", infile, N)

    except RuntimeError as e:
        print("Removing... ", infile)
        print(e)
        if not dry_run:
            os.remove(infile)

    except OSError as e:
        print("Removing... ",infile)
        print(e)
        if not dry_run:
            os.remove(infile)

    except KeyError as e:
        print("Removing... ",infile)
        print(e)
        if not dry_run:
            os.remove(infile)


#SBATCH --job-name=3D-targeted-extrusion-test
#SBATCH --metric=GPU
#SBATCH --time=05:00:00
#SBATCH --output=serial_test_%j.log
#SBATCH --array=0-26

echo "Running script as job $JOB_ID with task id: $TASK_ID on GPU: $GPU_INDEX"

DATAPATH="/net/levsha/share/agalicina/simulations/chromatin_fountains/"

FILELIST=($DATAPATH/data/traj1d/*.h5py)

INFILE="${FILELIST[$TASK_ID]}"
OUTFILE=${INFILE/traj1d/traj3d/}
OUTFILE=${OUTFILE/.h5py//}
/home/agalicina/anaconda3/bin/python -u $DATAPATH/01_3D_simulations.py -v -f -G $GPU_INDEX -i $INFILE -o $OUTFILE

#SBATCH --job-name=3D-targeted-extrusion-test
#SBATCH --metric=GPU
#SBATCH --time=20:00:00
#SBATCH --output=logs/serial_test_%j.log
#SBATCH --array=0-100

echo "Running script as job $JOB_ID with task id: $TASK_ID on GPU: $GPU_INDEX"

DATAPATH="/net/levsha/share/agalicina/simulations/chromatin_fountains/"

FILELIST=($(ls -t $DATAPATH/data/traj1d/*.h5py))
#($DATAPATH/data/traj1d/*.h5py)

INFILE="${FILELIST[$TASK_ID]}"
OUTFILE=${INFILE/traj1d/traj3d}
OUTFILE=${OUTFILE/.h5py/}

FILE="${OUTFILE}/blocks_10000-10099.h5"
if [ ! -f "$FILE" ]
then
	/home/agalicina/anaconda3/bin/python -u $DATAPATH/01_3D_simulations.py -v -f -G $GPU_INDEX -i $INFILE -o $OUTFILE
fi


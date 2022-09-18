#SBATCH --job-name=3D-targeted-extrusion-test
#SBATCH --metric=GPU
#SBATCH --time=05:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
#SBATCH --array=0

echo "Running script as job $JOB_ID with task id: $TASK_ID on GPU: $GPU_INDEX"

PATH="/net/levsha/share/agalicina/simulations/chromatin_fountains/"

/home/agalicina/anaconda3/bin/python -u python $PATH/01_3D_simulations.py -v -f -G $GPU_INDEX -i $PATH/data/traj1d/sea-rises_L100_S100_T50.h5py -o $PATH/data/traj3d/sea-rises_L100_S100_T50 

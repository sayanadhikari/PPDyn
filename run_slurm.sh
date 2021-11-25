#!/bin/bash

###############################################
# Script example for a normal MPI job on Fram #
###############################################

#SBATCH --job-name=ppdyndust01
#SBATCH --account=nn9299k
#SBATCH --time=7-00:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=4G

# set -o errexit # Exit the script on any error
# set -o nounset # Treat any unset variables as an error

module --quiet purge
module load Anaconda3/2020.11
export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null
conda activate /cluster/home/gauteah/.conda/envs/ppdyn
module list

export HDF5_USE_FILE_LOCKING='FALSE'

python ppdyn.py -i input.ini

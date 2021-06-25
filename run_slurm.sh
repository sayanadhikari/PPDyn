#!/bin/bash

###############################################
# Script example for a normal MPI job on Fram #
###############################################

#SBATCH --job-name=ppdyn01
#SBATCH --account=nn9299k
#SBATCH --time=0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30

set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variables as an error

module --quiet purge
module load Anaconda3/2019.07
module list

export HDF5_USE_FILE_LOCKING='FALSE'

# conda env create -f environment.yml
# conda activate ppdyn  

# pip install -r requirements.txt --user

srun python src/main.py -i input.ini

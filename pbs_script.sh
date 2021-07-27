#!/bin/bash

## JOB NAME
#PBS -N ppdyn_test
## QUEUE NAME
#PBS -q debugq
## COMPUTE RESOURCES REQUESTED FOR THE JOB SELECT = NO. OF CHUNKS/NODES, NCPUS = NO. OF CORES PER CHUNK/NODE
#PBS -l select=1:ncpus=40
## SPECIFY THE EXECUTION TIME LIMIT FOR THE CODE/APPLICATION IN HRS:MINS:SECS FORMAT
#PBS -l walltime=00:30:00
## JOIN THE OUTPUT AND ERROR FILES INTO A SINGLE FILE WITH NAME <JOBNAME>.O<JOBID>
#PBS -j oe
## EXPORT ALL ENVIRONMENT VARIABLES
#PBS -V
#EMAIL IS SENT WHEN THE JOB STARTS, TERMINATES AND ABORTS
#PBS -m bea
## SPECIFY EMAIL ADDRESS FOR NOTIFICATIONS
#PBS -M rinku.mishra@cppipr.res.in


# LOAD A MODULE BASED ON APPLICATION/CODE REQUIREMENT
module load anaconda/3
conda deactivate &>/dev/null
conda activate /home/rinku.mishra/.conda/envs/ppdyn
module list
export HDF5_USE_FILE_LOCKING='FALSE'

# WORKING DIRECTORY OF CODE/APPLICATION
#cd /scratch/scratch_run/deepakagg/queue_testing
# ENVIRONMENT VARIABLE FOR ACCESSING THE WORKING DIRECTORY WITH PBS VARIABLE
cd $PBS_O_WORKDIR

# WRITE IN NODES.TXT FILE, THE NODES ON WHICH THE RUN HAS BEEN LAUNCHED
#cat ${PBS_NODEFILE} > /scratch/scratch_data/deepakagg/queue_testing/nodes.txt

# RUN COMMAND BASED ON CODE/APPLICATION
time python ppdyn.py -i input.ini
#time mpirun -np 40 --machinefile $PBS_NODEFILE hostname > /scratch/scratch_data/deepakagg/queue_testing/output.txt

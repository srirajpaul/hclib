#!/bin/bash
##SBATCH -q scavenger
##SBATCH -q regular
#SBATCH -q debug
#SBATCH -N 2
##SBATCH -N 4
##SBATCH -N 8
##SBATCH -N 16
##SBATCH -N 32
##SBATCH -N 64
#SBATCH -C haswell
##SBATCH -t 01:59:00
#SBATCH -t 00:30:00
##SBATCH --time-min=03:30:00

#N=2
#N=4
#N=8
#N=16
#N=32
#N=64

## This script will run the experiment for the number of nodes requested using 'SBATCH -N'
## To select the jos queue "SBATCH -q', number of nodes 'SBATCH -N' and time 'SBATCH -t'
## Repeat the experiment by setting 'SBATCH -N 2' then 4, 8, 16, 32, 64

date

N=$SLURM_JOB_NUM_NODES

/bin/bash run.sh -s l -n 5 -r 32 -t $N -m $N -u 1 |& tee out${N}.log

date


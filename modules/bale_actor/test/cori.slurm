#!/bin/bash
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

#N=$SLURM_JOB_NUM_NODES
N=64

srun -n $N ./ig_conveyor 
srun -n $N ./ig_selector


srun -n $N ./ig_conveyor -n 10000000
srun -n $N ./ig_selector -n 10000000


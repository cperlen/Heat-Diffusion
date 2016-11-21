#!/bin/bash
# Parallel job using 16 processors:
#SBATCH -N 1 # number of nodes
#SBATCH --ntasks-per-node=16 # number of processors per node
#SBATCH -t 0:15:00 # run for 15 minutes max
# Load openmpi environment
module load openmpi


for nx in 128 256 512
do
srun ./heat_mpi $nx > heat_mpi.$nx.16.out
done
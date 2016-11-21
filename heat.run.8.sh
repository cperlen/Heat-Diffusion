#!/bin/bash
# Parallel job using 8 processors:
#SBATCH -N 1 # number of nodes
#SBATCH --ntasks-per-node=8 # number of processors per node
#SBATCH -t 0:15:00 # run for 15 minutes max
# Load openmpi environment
module load openmpi


for nx in 128 256 512
do
./heat_omp $nx 8 > heat_omp.$nx.8.out
srun ./heat_mpi $nx > heat_mpi.$nx.8.out
done
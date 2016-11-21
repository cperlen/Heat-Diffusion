#!/bin/bash
# Parallel job using 4 processors:
#SBATCH -N 1 # number of nodes
#SBATCH --ntasks-per-node=4 # number of processors per node
#SBATCH -t 0:15:00 # run for 15 minutes max
# Load openmpi environment
module load openmpi


for nx in 128 256 512
do
./heat_omp $nx 4 > heat_omp.$nx.4.out
srun ./heat_mpi $nx > heat_mpi.$nx.4.out
done
#!/bin/bash -l

#BATCH -A g2017012

#SBATCH -t 8:00

#SBATCH -p node -n 36


module load gcc openmpi
make

mpirun -np 1 ./MPI_wave 200 1 1

mpirun -np 4 ./MPI_wave 200 2 2

mpirun -np 9 ./MPI_wave 200 3 3

mpirun -np 16 ./MPI_wave 200 4 4

mpirun -np 25 ./MPI_wave 200 5 5

mpirun -np 36 ./MPI_wave 200 6 6

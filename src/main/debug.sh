#!/bin/bash
#SBATCH --partition=debug                # Name of Partition
#SBATCH --ntasks=8                                # Request of CPU cores
#SBATCH --time=00:30:00                              # Job should run for up to 1.5 hours (for example)
#SBATCH --mail-type=END                              # Event(s) that triggers email notification (BEGIN,END,FAIL,ALL)
##SBATCH --mail-user=<huilin.ye@uconn.edu>           # Destination email address

srun --mpi=openmpi palammps in.adsphere >log.file


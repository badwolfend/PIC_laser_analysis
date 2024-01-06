#!/bin/bash
#SBATCH -J /scratch/jyoung36/osiris-1.0.0/bin/osiris--2D.e
#SBATCH -p preempt 
#SBATCH -o job_osiris_output_%j
#SBATCH --mail-type=all
#SBATCH --mem-per-cpu=1GB
#SBATCH -t 1:59:00
#SBATCH -n 4 
#SBATCH -N 1
#SBATCH -w bhd0060
hostname
module load intel impi hdf5 gcc git
mpiexec -np $SLURM_NTASKS ../../bin/osiris--2D.e laser-2d


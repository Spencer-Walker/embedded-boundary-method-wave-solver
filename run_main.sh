#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=00:30:00
#SBATCH --partition=shas-testing
#SBATCH --job-name=eb-pml-job
#SBATCH --output=eb-pml-job-%j.out
#SBATCH --mail-type=start,end
#SBATCH --mail-user=peter.rosenthal@colorado.edu

module purge
module load intel impi petsc python matlab

python run.py

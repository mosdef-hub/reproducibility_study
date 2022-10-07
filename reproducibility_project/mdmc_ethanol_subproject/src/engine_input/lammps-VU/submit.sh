#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=96:00:00
#SBATCH --partition=week-long-std
#SBATCH --output=test_job_out.txt
#SBATCH --export=ALL

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
echo $infile $seed $T $P $rcut $tstep
module load lammps/3Aug2022
module avail
lmp -in $infile -var seed $seed -var T $T -var P $P -var rcut $rcut -var tstep $tstep

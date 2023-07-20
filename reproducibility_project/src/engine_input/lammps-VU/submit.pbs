#!/bin/bash -l
#PBS -N FILL_IN_NAME
#PBS -m ae
#PBS -M nicholas.craven.76@gmail.com
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -q standard
#PBS -j oe
cd $PBS_O_WORKDIR
module load lammps
export OMP_NUM_THREADS=16
mpirun -np 8 lmp < $infile -var seed $seed -var T $T -var P $P -var rcut $rcut -var tstep $tstep
echo $PBS_JOBID $PBS_O_WORKDIR >> ~/job-id-dirs.txt

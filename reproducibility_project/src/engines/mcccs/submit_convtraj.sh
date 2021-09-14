#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 200:00:00
#SBATCH --mem=8G
#SBATCH -o outconvtraj.out
#SBATCH -e errconvtraj.err
#SBATCH --job-name=conv_traj

cd $SLURM_SUBMIT_DIR
conda --version
#source activate /home/siepmann/singh891/.conda/envs/halogen37
conda activate mosdef-study38
date
python conv_traj.py
date

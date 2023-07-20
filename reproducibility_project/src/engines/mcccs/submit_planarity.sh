#!/bin/sh
#SBATCH --job-name=planarity
#SBATCH -t 199:59:59
#SBATCH --ntasks=1
#SBATCH --mem=8g
module purge
conda --version
source /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate mosdef-study38
python Planarity.py

#!/bin/sh
#SBATCH --job-name=bl_analysis
#SBATCH -t 199:59:59
#SBATCH --ntasks=1
#SBATCH --mem=4g
module purge
conda --version
source /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate mosdef-study38
python bl_analysis_overall.py

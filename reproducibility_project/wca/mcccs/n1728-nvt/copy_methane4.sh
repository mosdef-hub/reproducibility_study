#!/bin/bash
#SBATCH --job-name=copy_methane4    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request

#SBATCH --time=24:05:00               # Time limit hrs:min:sec

#SBATCH --output=parallel_%j.log     # Standard output and error log

. /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate mosdef-study38

rsync -av /home/rs/space/projects/final_repro_methane/reproducibility_study/reproducibility_project/methane_systemsize_subproject4/* .

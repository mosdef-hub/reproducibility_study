{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

#SBATCH --nodes=1
#SBATCH --partition=week-long-tesla
#SBATCH --time=72:00:00
#SBATCH --gres gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=co.d.quach@vanderbilt.edu
#SBATCH -o output-%j.log
#SBATCH -e error-%j.log

module load anaconda/3.9
module load hoomd/single/4.0.0tesla

export MOSDEF_PYTHON=~/.conda/envs/hoomd4.0.0quachcd/bin/python
conda activate --stack mosdef-study38

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}

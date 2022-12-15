{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% if gpus %}
#SBATCH -q standard
#SBATCH --gres gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=co.d.quach@vanderbilt.edu
#SBATCH -o output-%j.log
#SBATCH -e error-%j.log

conda activate mosdef-study38

module load gromacs/2020.6

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}

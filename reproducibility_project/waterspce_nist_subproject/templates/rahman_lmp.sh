{% extends "slurm.sh" %}
{% block header %}
#!/bin/sh
#SBATCH --job-name="{{ id }}"
#SBATCH -t 24:00:00
#SBATCH --partition=day-long-std
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=2g
#SBATCH --output=test_job_out.txt
#SBATCH --export=ALL

module load lammps/3Aug2022
{% endblock header %}

{% block body %}
	{{- super () -}}

{% endblock body %}

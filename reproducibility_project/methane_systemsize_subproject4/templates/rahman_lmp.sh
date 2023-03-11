{% extends "slurm.sh" %}
{% block header %}
#!/bin/sh
#SBATCH --job-name="{{ id }}"
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2g
{% endblock %}

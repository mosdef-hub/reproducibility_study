{% extends "slurm.sh" %}
{% block header %}
#!/bin/sh
#SBATCH --job-name="{{ id }}"
#SBATCH -t 199:59:59
#SBATCH --ntasks=1
#SBATCH --mem=1G
module purge
module load anaconda
#conda --version
#source /home/rs/anaconda3/etc/profile.d/conda.sh
source activate mosdef-study38
date >> execution.log

{% if partition %}
#SBATCH --partition={{ partition }}
{% endif %}
{% if nodelist %}
#SBATCH --nodelist={{ nodelist }}
{% endif %}
#SBATCH --output=workspace/{{operations[0]._jobs[0]}}/job_%j.o
#SBATCH --error=workspace/{{operations[0]._jobs[0]}}/job_%j.e
{% endblock %}

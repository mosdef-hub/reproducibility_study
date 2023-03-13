{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q standard
#PBS -m abe
#PBS -M quachcd.rahman@gmail.com

module load anaconda/v3.6_5.2
source activate mosdef-study38
module load gromacs/2020.6

{% endblock header %}

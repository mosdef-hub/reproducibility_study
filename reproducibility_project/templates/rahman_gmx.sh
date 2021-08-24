{% extends base_script %}
{% block project_header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes={{ nn }}:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q low
#PBS -m abe
#PBS -M quachcd.rahman@gmail.com

module load anaconda/v3.6_5.2
source activate mosdef-study38

{% endblock project_header %}

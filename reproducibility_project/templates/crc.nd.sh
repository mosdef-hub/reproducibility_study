{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash
#$ -N {{ id }}
#$ -pe smp {{ np_global }}
#$ -r n
#$ -q long
#$ -m ae
#$ -M rsmith56@nd.edu

conda activate mosdef-study38
module load ompi
export PATH=/afs/crc.nd.edu/user/r/rsmith56/software/Cassandra_master:${PATH}

{% block tasks %}
{% endblock %}
{% endblock %}

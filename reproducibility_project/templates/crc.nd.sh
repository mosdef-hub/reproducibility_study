{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash
#$ -N {{ id }}
#$ -pe smp {{ np_global }}
#$ -r n
#$ -q long
#$ -m abe
#$ -M emarinri@nd.edu

conda activate mosdef-study38

{% block tasks %}
{% endblock %}
{% endblock %}

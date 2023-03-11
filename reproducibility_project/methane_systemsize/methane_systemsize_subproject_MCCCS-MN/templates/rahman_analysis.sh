{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q standard
#PBS -m abe
#PBS -M gilmerjb.job.scheduler@gmail.com

conda activate mosdef-study38

{% endblock header %}

{% block body %}
{% set cmd_suffix = cmd_suffix|default('') ~ (' &' if parallel else '') %}
{% for ops_batch in operations|batch(16) %}
{% for operation in ops_batch %}
{{ operation.cmd }} || : {{ cmd_suffix }}
{% endfor %}
wait
{% endfor %}

{% endblock body %}

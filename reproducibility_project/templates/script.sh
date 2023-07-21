{% extends "delta.sh" %}

{% block header %}
    {{- super () -}}
{% endblock header %}
{% block custom_content %}
{#
    This block is not used by any other template and can be safely modified
    without the need to call super(). We recommend most additions to the
    templates go here if they are not direct changes to an existing template.

    For example, commands like `module load ...` or printing diagnostic
    information from the scheduler can be done in this block.
#}
source $HOME/mamba-env.sh
conda activate mosdef-study38
export MOSDEF_PYTHON=python

{% endblock custom_content %}
{% block body %}
    {{- super () -}}
{% endblock body %}

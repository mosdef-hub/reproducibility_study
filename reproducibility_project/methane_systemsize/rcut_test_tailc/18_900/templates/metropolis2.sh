#SBATCH -t 199:59:59
module purge
conda --version
source /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate mosdef-study38
date >> execution.log
{{ super() }}
{% endblock %}

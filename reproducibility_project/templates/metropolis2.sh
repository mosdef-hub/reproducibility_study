module purge
#SBATCH -t 200:00:00
conda --version
source /home/rs/anaconda3/etc/profile.d/conda.sh
conda activate mosdef-study38
date >> execution.log
{{ super() }}
{% endblock %}

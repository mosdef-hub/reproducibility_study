{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

#SBATCH --nodes=1
#SBATCH --partition=week-long-tesla
#SBATCH --time=72:00:00
#SBATCH --gres gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=craven76@gmail.com
#SBATCH -o output-%j.log
#SBATCH -e error-%j.log

echo "Running hoomd"
module load anaconda/3.9
#module load hoomd/single/3.11.0tesla
module load hoomd/single/4.0.0tesla
echo $(conda list hoomd)
echo $(which python)
echo "Printing HOOMD"
#module load hoomd/single/4.0.0tesla
#module load hoomd/single/2.9.7tesla

#export MOSDEF_PYTHON=~/.conda/envs/hoomd3.11.0.quachcd/bin/python
#conda activate hoomd3.11.0.quachcd
conda activate hoomd3.11.0.craven76
export MOSDEF_PYTHON=~/.conda/envs/hoomd3.11.0.craven76/bin/python
echo $(conda list hoomd)
echo $(which python)
echo $(conda info)
#export MOSDEF_PYTHON=~/.conda/envs/craven76.hoomd4/bin/python
#conda activate craven76.hoomd4

#export MOSDEF_PYTHON=~/anaconda3/envs/hoomd4.0.0.craven76/bin/python
#echo $MOSDEF_PYTHON
#echo $(conda info --envs)
#conda activate /raid6/homes/craven76/anaconda3/envs/hoomd4.0.0.craven76


{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}

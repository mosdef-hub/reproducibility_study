# Reproducibility Study
Repo for data collection, discussion, etc for a [MoSDeF](http://mosdef.org) reproducibility study.

## Installation
We use [conda](https://docs.conda.io/en/latest/) to manage our software environment, so you'll need a conda installation--we recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html) with [mamba](https://github.com/mamba-org/mamba).

Once conda and mamba are installed, the `mosdef-study38` environment can be created by running:

```bash
mamba env create -f environment.yml
```
And the environment can be activated by running:
```bash
conda activate mosdef-study38
```

### Installing simulation engines
The simulation engines used in this study require some additional setup.
GROMACS, LAMMPS, HOOMD-blue, and Cassandra can be add to the environment by running:
```bash
mamba env update -f engines.yml
```
NOTE: For hoomd this installation may not work unless executed on a node with a GPU. To check GPU configuration run `nvidia-smi`. The hoomd gpu api version (`python -c "import hoomd; print(hoomd.version.gpu_api_version)"`) should be the same as the cuda version.

GOMC can be installed according to its documentation [here](https://gomc.eng.wayne.edu/).

MCCCS can be installed according to its documentation [here](https://ccs-psi.org/node/52).

## Use
This project uses the [Signac framework](https://signac.io/) to manage its parameter space. Instructions for initializing and submitting/running the simulation and analysis workflows can be found in [the project guide](reproducibility_project/README.md).

## Sharing data
In order to share data between contributors, we use [rclone](https://rclone.org/) to sync to a shared [dropbox](https://www.dropbox.com) folder.

rclone can be installed without sudo permissions (as might be necessary, e.g., on a computing cluster) using conda:
```
conda install -c conda-forge rclone
```

Next [rclone must be linked to dropbox](https://rclone.org/dropbox/). The next steps will use the remote name specified during the linking process. For the following examples we will use the name `dropbox`. If you forgot the name, you can use `rclone listremotes`.

Before syncing the workspace data, it is recommended to remove any empty job folders, as these will overwrite the data on dropbox. One way to do this is to remove the job folders of unused engines. We provide a bash script which can help(located in `reproducibility_study/reproducibility_project/bin/clean_by_engine.sh`). For example, to remove job folders for all engines EXCEPT hoomd, run:

```
bash bin/clean_by_engine.sh hoomd
```

Be careful when using this script as it could delete your data! Make sure so check which project (`mosdef_reproducibility` or `lrc_shift_subproject`) you are in using `signac project`. Before using it is advisable to make a backup and make sure that you spell/capitalize your engine the same way as in the statepoint.

Once your workspace contains only the data you want to sync, you can see how rclone would copy it to dropbox using the `--dry-run` flag:

```
rclone copy workspace dropbox:MoSDeF\ Repro\ Study\ Data/workspace --dry-run
```

The same command without `--dry-run` can be used to copy the data over.

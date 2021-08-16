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

GOMC can be installed according to its documentation [here](https://gomc.eng.wayne.edu/).

MCCCS can be installed according to its documentation [here](https://ccs-psi.org/node/52).

## Use
This project uses the [Signac framework](https://signac.io/) to manage its parameter space. Instructions for initializing and submitting/running the simulation and analysis workflows can be found in [the project guide](project/README.md).

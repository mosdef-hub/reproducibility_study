# Reproducibility Study
Repo for data collection, discussion, etc for a [MoSDeF](http://mosdef.org) reproducibility study.

## Installation
We use [conda](https://docs.conda.io/en/latest/) to manage our software environment, so you'll need a conda installation--we recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html) with [mamba](https://github.com/mamba-org/mamba).

Once conda and mamba are installed, the `mosdef-study38` environment can be created by running:

```bash
mamba env create -f environment.yml
```

Each simulation engine requires additional setup.
TODO engine-specific install steps

## Use
This project uses the [Signac framework](https://signac.io/) to manage its parameter space. Instructions for initializing and submitting/running the simulation and analysis workflows can be found in [the project guide](project/README.md).

### Creating useful git commit messages
For most additions to this project, we expect users to provide useful/detailed commit messages.
Please refer to the commit messages in the current `git log` for examples.
[This is a great resource for creating detailed commit messages.](https://chris.beams.io/posts/git-commit/)


### Using `pre-commit` to auto-format commits

**NOTE** We cannot set up `pre-commit-ci` to auto format this repository at the moment since it is a private repo.
`pre-commit-ci` is free for **public** repositories, but not private ones like this project.

Users are expected to set up `pre-commit` on their local machines for the time being.

```bash
conda activate mosdef-study38
mamba install -c conda-forge pre-commit
pre-commit install --install-hooks
pre-commit run --all-files
```

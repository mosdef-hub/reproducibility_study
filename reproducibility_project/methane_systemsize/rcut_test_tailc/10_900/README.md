## You are here
Methane systemsize subproject for reproducibility study.

## Objective
Comparing density results from various system sizes to see the effect of system size on densit values.

## How to set up the simulations
This subproject can be run from the `subproject-init.py`. Once executed, the jobs for methane NpT simulations will be created for two engines, lammps-VU and mcccs. These engines will demonstrate the effect of system size on density values obtained from MD and MC engines. The densities for three system sizes are compared between engines to investigate the effect of COM degrees of freedom removal in MD on system density.



## How to generate this project data from scratch

**IMPORTANT: All the following commands are intended to be run from this location with the `mosdef-study38` environment active.**

To initialize the project workspace:
```bash
python subproject-init.py
```
You will see the empty workspace folder populated with the hashed job directories.

Each engine directory in `src/engines` contains its own `project.py` file. To query the status of simulations for a specific engine (in the following example we use hoomd), you would run:
```bash
python src/engines/hoomd/project.py status
```

To learn more about the ways to use this project:
```bash
python src/engines/hoomd/project.py -h
```
or check out the [FlowProject documentation](https://docs.signac.io/en/latest/flow-project.html).


## Layout

```
├── README.md
├── __init__.py
├── subproject-init.py
├── signac.rc
├── signac_project_document.json
├── templates
├── analysis-ethanol-bonds.py
├── src
│   ├── __init__.py
│   ├── engine_input
│   │   ├── __init__.py
│   │   ├── lammps
│   │   └── mcccs
│   ├── engines
│   │   ├── __init__.py
│   │   ├── lammps
│   │   │   └── project.py
│   │   └── mcccs
│   │       └── project.py
└── workspace
```

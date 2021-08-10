## You are here
Main location for reproducibility study.

## How to use
All commands are intended to be run from this location with the `mosdef-study38` environment active.
To initialize the project workspace:
```bash
python init.py
```
You will see the empty workspace folder populated with the hashed job directories.

Each engine directory in `src/engines` contains its own `project.py` file. To run the simulations for a specific engine (in the following example we use hoomd), you would run:
```bash
python src/engines/hoomd/project.py run
```

Finally, to run the analysis of the simulation output you would run:
```bash
python project-analysis.py run
```
## Dashboard instructions
(TODO)

## Layout

```
├── README.md
├── __init__.py
├── init.py
├── project-analysis.py
├── signac.rc
├── signac_project_document.json
├── src
│   ├── README.md
│   ├── __init__.py
│   ├── dashboard.py
│   ├── engine_input
│   │   ├── __init__.py
│   │   ├── cassandra
│   │   ├── gomc
│   │   ├── gromacs
│   │   │   └── mdp
│   │   │       ├── em.mdp
│   │   │       ├── npt.mdp
│   │   │       └── nvt.mdp
│   │   ├── hoomd
│   │   ├── lammps
│   │   └── mcccs
│   │       ├── fort77maker_onebox.py
│   │       ├── methane
│   │       │   ├── fort.4.cool
│   │       │   ├── fort.4.equil
│   │       │   ├── fort.4.melt
│   │       │   ├── fort.4.prod
│   │       │   └── topmon.inp
│   │       └── pentane
│   │           ├── fort.4
│   │           └── topmon.inp
│   ├── engines
│   │   ├── README.md
│   │   ├── __init__.py
│   │   ├── cassandra
│   │   │   └── project.py
│   │   ├── gomc
│   │   │   └── project.py
│   │   ├── gromacs
│   │   │   └── project.py
│   │   ├── hoomd
│   │   │   ├── project.py
│   │   ├── lammps
│   │   │   └── project.py
│   │   └── mcccs
│   │       ├── fort77maker_onebox.py
│   │       └── project.py
│   ├── molecules
│   │   ├── __init__.py
│   │   └── methane_ua.py
│   └── xmls
│       └── spce.xml
├── tests
│   ├── __init__.py
│   ├── base_test.py
│   ├── test_methane_ua.py
│   └── test_spce_water.py
└── workspace
```

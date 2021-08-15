## You are here
Main location for reproducibility study.

## How to use

**IMPORTANT: All the following commands are intended to be run from this location with the `mosdef-study38` environment active.**

To initialize the project workspace:
```bash
python init.py
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

Finally after running or submitting the engine-specific simulation projects, to run the analysis of the simulation output you would run:
```bash
python project-analysis.py run
```

## Dashboard instructions
[Signac-dashboard](https://docs.signac.io/projects/dashboard/en/latest/) is a convenient application for displaying a signac project.

To initialize the dashboard:
```bash
python src/dashboard.py run
```

This will start the flask app and output the following message to the screen:
```
 * Serving Flask app 'signac-dashboard' (lazy loading)
 * Environment: production
   WARNING: This is a development server. Do not use it in a production deployment.
   Use a production WSGI server instead.
 * Debug mode: off
 * Running on http://localhost:8888/ (Press CTRL+C to quit)
```
Copy and paste the server and port address (`http://localhost:8888/`) into your browser to view the dashboard.

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

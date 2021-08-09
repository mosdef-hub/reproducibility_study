Main location for reproducibility study.

Layout is as below (subject to change):

```
├── README.md
├── __init__.py
├── __pycache__
│   └── __init__.cpython-37.pyc
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
│   ├── __pycache__
│   │   ├── __init__.cpython-37.pyc
│   │   ├── base_test.cpython-37-pytest-6.2.4.pyc
│   │   └── test_methane_ua.cpython-37-pytest-6.2.4.pyc
│   ├── base_test.py
│   ├── test_methane_ua.py
│   └── test_spce_water.py
└── workspace
```

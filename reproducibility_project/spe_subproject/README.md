## You are here
Single Point Energy (SPE) subproject for reproducibility study

## How to set up single energy calculations
Single point energies are calculated from the snapshots found in `src/system_snapshots`. These are created from the file `src/create_snapshots.py`. This utilizes the system builder in the main study to create initial configurations that can be compared across simulation engine for energy calculations.

**NOTE: Packmol is used to pack the simulations box, and can generate different packings with the same seed across Linux and macOS operating distributions. The simulation snapshots used in this study were created using macOS Catalina Version 10.15.7.**

The energy calculations will all be output to `log-spe.txt` and located in the job workspace directories. All energy units will be reported in **extensive**(no division by number of moleules) energies of kJ/mol.

The energy will be logged as space separated columns with the first row containing headers. They will be as follows:
total_energy potential_energy kinetic_energy vdw_energy coul_energy pair_energy bonds_energy angles_energy dihedrals_energy

The following energy definitions are calculated using the [lammps thermo](https://docs.lammps.org/thermo_style.html) calculation methods.
total_energy = etotal
potential_energy = pe
kinetic_energy = ke
vdw_energy = evdwl
coul_energy = ecoul
pair_energy = epair
bonds_energy = ebond
angles_energy = eangle
dihedrals_energy = edihed
**NOTE: Please log a blank column of data for any energies that are not output or calculated by your simulation engine**
This can be done by using the `None` value for the entries in that pandas dataframe. To write out your data, use:
```python
df.to_csv('log-spe.txt', header=True, index=False, sep=" ")
```

## How to use from scratch

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
├── src
│   ├── README.md
│   ├── __init__.py
│   ├── dashboard.py
│   ├── engine_input
│   │   ├── __init__.py
│   │   ├── cassandra
│   │   ├── gomc
│   │   ├── gromacs
│   │   ├── hoomd
│   │   ├── lammps
│   │   └── mcccs
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
│   │       └── project.py
|   ├── system_snapshots
│   │   ├── waterSPCE.json
│   │   ├── ethanolAA.json
│   │   ├── benzeneUA.json
│   │   ├── pentaneUA.json
│   │   ├── methaneUA.json
│   │   ├── create_snapshots.py
└── workspace
```

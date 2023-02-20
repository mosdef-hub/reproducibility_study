## You are here
Ethanol Bondlength subproject for reproducibility study.

## How to set up single energy calculations
This subproject can be run from the `subproject-init.py`. Once executed, the jobs will be created for two engines, lammps-VU and mcccs. These engines will demonstrate the effects of MD and MC engines with respective rigid and flexible OH bonds in the ethanol molecule. The resulting densities and average bond lengths are compared to the main project treatment of this molecule to demonstrate the systematic error seen between the MD and MC engines.

The energy will be logged as "comma" separated columns with the first row containing headers. They will be as follows:
potential_energy, tot_vdw_energy, tail_energy, tot_electrostatics, short_range_electrostatics, long_range_electrostatics, tot_pair_energy, bonds_energy, angles_energy, dihedrals_energy, tot_bonded_energy,
intramolecular_energy, intermolecular_energy.

## Reported Energy Definitions
The following energy definitions are calculated using the [lammps thermo](https://docs.lammps.org/thermo_style.html) calculation methods. </br>
```
* potential_energy = pe </br>
* tot_vdw_energy = evdwl </br>
* tail_energy = etail </br>
* tot_electrostatics = elong </br>
* short_range_electrostatics = ecoul </br>
* long_range_electrostatics = elong </br>
* tot_pair_energy = epair </br>
* bonds_energy = ebond </br>
* angles_energy = eangle </br>
* dihedrals_energy = edihed </br>
* tot_bonded_energy = emol </br>
* intramolecular_energy = emol + intramolecular ecoul + intramolecular evdw </br>
* intermolecular_energy = intermolecular ecoul + interolecular evdw + intermolecular elong </br>
```
</br>
**NOTE: Please log a blank column of data for any energies that are not output or calculated by your simulation engine**
This can be done by using the `None` value for the entries in that pandas dataframe. To write out your data, use:
```python
df.to_csv('log-spe.txt', header=True, index=False, sep=",")
```

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

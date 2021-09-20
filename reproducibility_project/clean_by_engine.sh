#!/usr/bin/env bash
echo Removing jobs from all engines but $1
for ENGINE in cassandra mcccs gomc gromacs hoomd lammps-VU lammps-UD
do
    if [[ $ENGINE != $1 ]]
    then
        signac rm $(signac find engine $ENGINE)
        echo Removed jobs for $ENGINE
    fi
done

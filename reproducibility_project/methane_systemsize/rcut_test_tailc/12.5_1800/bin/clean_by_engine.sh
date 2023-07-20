#!/usr/bin/env bash
echo Removing jobs from all engines but $1
for ENGINE in cassandra mcccs gomc gromacs hoomd lammps-VU lammps-UD
do
    if [[ $ENGINE != $1 ]]
    then
        read -p "Delete $ENGINE job directories?[y/N] " yn
        case $yn in
            [Yy]* ) signac rm $(signac find engine $ENGINE); echo Removed jobs for $ENGINE;;
            * ) echo Skipping $ENGINE;;
        esac
    fi
done

"""Move and remove obsolete jobs from workspace. Use with care."""
import json
import os
import os.path
import shutil

import signac
import unyt as u


def clean_workspace():
    """Move and remove obsolete jobs from workspace. Use with care."""
    if not os.path.isdir("archives"):
        os.mkdir("archives")

    project = signac.get_project()
    if project.num_jobs() == 992:
        raise Exception("Number of jobs is 992 (expected). Abort cleaning.")

    if not os.path.isdir("archives/lammps-UD_workspace"):
        os.mkdir("archives/lammps-UD_workspace")

    lammpsUD_jobs = list()
    for job in project.find_jobs({"engine": "lammps-UD"}):
        lammpsUD_jobs.append(job.id)

    for job in lammpsUD_jobs:
        try:
            shutil.move(
                f"workspace/{job}", f"archives/lammps-UD_workspace/{job}"
            )
        except:
            if os.path.isdir(f"archives/lammps-UD_worspace"):
                print("workspace/{job} has already been moved")
            else:
                print(f"workspace/{job} not found.")

    # Remove obsolete molecule nomenclature
    standard_molecules = [
        "methaneUA",
        "pentaneUA-flexible_bonds",
        "pentaneUA-constrain_bonds",
        "benzeneUA",
        "waterSPCE",
        "ethanolAA",
    ]

    for job in project:
        if job.sp.molecule not in standard_molecules:
            job.remove()

    # Move duplicated MCCCS jobs
    updated_masses = {
        "methaneUA": 16.043,
        "pentaneUA-flexible_bonds": 72.151,
        "pentaneUA-constrain_bonds": 72.151,
        "benzeneUA": 78.114,
        "waterSPCE": 18.015,  # 18.015324,
        "ethanolAA": 46.069,  # 46.068672,
    }
    obsolete_mcccs = list()
    for job in project:
        if (
            job.sp.engine == "mcccs"
            and job.sp.mass != updated_masses[job.sp.molecule]
        ):
            obsolete_mcccs.append(job.id)

    # Move obsolete MCCCS jobs
    if not os.path.isdir("archives/mcccs_workspace"):
        os.mkdir("archives/mcccs_workspace")

    for mcccs_job in obsolete_mcccs:
        try:
            shutil.move(
                f"workspace/{mcccs_job}", f"archives/mccc_workspace/{mcccs_job}"
            )
        except:
            if os.path.isdir(f"archives/mcccs_workspace/{mcccs_job}"):
                print(f"workspace/{mcccs_job} has already been moved.")
            else:
                print(f"workspace/{mcccs_job} not found.")

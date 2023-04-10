"""Copy LAMMPS-VU data from the main aggregate."""

import os
import shutil

import signac


def copy_lammps_etoh_flex():
    """Copy the lammps job with flexible OH bond."""
    main_aggregate_project = signac.get_project("../../aggregate_summary/")
    assert main_aggregate_project.id == "aggregate_main_summary"

    mdmc_aggregate_project = signac.get_project()
    assert mdmc_aggregate_project.id == "aggregate_mdmc_etoh_summary"

    for job in main_aggregate_project.find_jobs(
        {"engine": "lammps-VU", "molecule": "ethanolAA"}
    ):
        shutil.copytree(
            job.ws, f"{mdmc_aggregate_project.workspace()}/{job.id}"
        )

    job_ids = set()
    for job in mdmc_aggregate_project.find_jobs(
        {"engine": "lammps-VU", "ensemble": "NPT"}
    ):
        job_ids.add(job.id)
        job.update_statepoint({"ensemble": "NPT-flexOH"}, overwrite=True)

    new_job_ids = set()
    for job in mdmc_aggregate_project.find_jobs(
        {"engine": "lammps-VU", "ensemble": "NPT-flexOH"}
    ):
        new_job_ids.add(job.id)

    assert job_ids != new_job_ids
    assert len(job_ids) == len(new_job_ids)

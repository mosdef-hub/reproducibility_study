"""Setup for signac, signac-flow, signac-dashboard for this study."""
import pathlib
from typing import List

import flow
import numpy as np
import signac
from flow import aggregator


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


def _determine_sampling_information(
    job: signac.contrib.project.Job,
    ensemble: str,
    prop: str,
    filename: str = None,
) -> None:
    """Write out sampling results for production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    if filename is None:
        filename = f"log-{ensemble}.txt"
    sample_job(
        job,
        ensemble=ensemble,
        filename=filename,
        variable=prop,
        threshold_fraction=0.75,
        threshold_neff=100,
    )


def _is_prop_subsampled(
    job: signac.contrib.Project.Job, ensemble: str, prop: str
):
    """Check if the property has been subsampled."""
    return job.doc.get(f"{ensemble}/sampling_results", {}).get(f"{prop}", None)


def _get_largest_t0(
    job: signac.contrib.Project.Job, ensemble: str, props: List[str]
) -> int:
    """Return the largest t0 value for all sampling results."""
    index_list = list()
    for prop in props:
        prop_str = f"{ensemble}/sampling_results"
        prop_dict = job.doc.get(prop_str, {}).get(prop, None)
        if prop_dict is not None:
            index_list.append(prop_dict["start"])
    return max(index_list)


def _calc_statistics(
    job: signac.contrib.Project.Job, ensemble: str, prop: str
) -> None:
    """Calculate avg and std of subsampled data."""
    with job:
        with job.stores[f"{ensemble}_{prop}"] as data:
            job.doc[f"{ensemble}_{prop}_avg"] = np.mean(data["property"])
            job.doc[f"{ensemble}_{prop}_std"] = np.std(data["property"])


def all_npt_props_averaged(job: signac.contrib.Project.Job):
    """Check if all npt properties are averaged."""
    props = [
        "potential_energy",
        "kinetic_energy",
        "temperature",
        "pressure",
        "density",
    ]
    truthy = list()
    for prop in props:
        prop_string = f"npt_{prop}"
        is_value = [
            job.doc.get(prop_string + "_avg", False),
            job.doc.get(prop_string + "_std", False),
        ]
        truthy.append(all(is_value))
    return all(truthy)


@Project.label
def npt_ke_subsampled(job):
    """Check if kinetic energy has been subsampled."""
    return _is_prop_subsampled(job, ensemble="npt", prop="kinetic_energy")


@Project.label
def nvt_ke_subsampled(job):
    """Check if kinetic energy has been subsampled."""
    return _is_prop_subsampled(job, ensemble="nvt", prop="kinetic_energy")


@Project.label
def npt_pe_subsampled(job):
    """Check if potential energy has been subsampled."""
    return _is_prop_subsampled(job, ensemble="npt", prop="potential_energy")


@Project.label
def nvt_pe_subsampled(job):
    """Check if potential energy has been subsampled."""
    return _is_prop_subsampled(job, ensemble="nvt", prop="potential_energy")


@Project.label
def npt_density_subsampled(job):
    """Check if density has been subsampled."""
    return _is_prop_subsampled(job, ensemble="npt", prop="density")


@Project.label
def npt_pressure_subsampled(job):
    """Check if pressure has been subsampled."""
    return _is_prop_subsampled(job, ensemble="npt", prop="pressure")


@Project.label
def npt_temp_subsampled(job):
    """Check if temperature has been subsampled."""
    return _is_prop_subsampled(job, ensemble="npt", prop="temperature")


@Project.label
def nvt_temp_subsampled(job):
    """Check if temperature has been subsampled."""
    return _is_prop_subsampled(job, ensemble="nvt", prop="temperature")


@Project.label
def nvt_volume_subsampled(job):
    """Check if volume has been subsampled."""
    return _is_prop_subsampled(job, ensemble="nvt", prop="volume")


@Project.label
def npt_prod_finished(job):
    """Generate label if npt production is complete."""
    return job.isfile("trajectory-npt.gsd") and job.isfile("log-npt.txt")


@Project.label
def nvt_prod_finished(job):
    """Generate label if nvt production is complete."""
    return job.isfile("trajectory-nvt.gsd") and job.isfile("log-nvt.txt")


@Project.label
def npt_rdf_calculated(job):
    """Generate label if npt rdf calc is complete."""
    return (
        job.isfile("trajectory-npt.gsd")
        and job.isfile("npt_rdf.png")
        and job.isfile("npt_rdf.txt")
    )


@Project.label
def nvt_rdf_calculated(job):
    """Generate label if nvt rdf calc is complete."""
    return (
        job.isfile("trajectory-nvt.gsd")
        and job.isfile("nvt_rdf.png")
        and job.isfile("nvt_rdf.txt")
    )


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
@Project.post(lambda job: job.isfile("npt_rdf.png"))
@Project.post(lambda job: job.isfile("npt_rdf.txt"))
def rdf_npt_analysis(job):
    """Run RDF analysis."""
    from reproducibility_project.src.analysis.rdf import gsd_rdf

    # RDF
    gsd_rdf(job, ensemble="npt")


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
@Project.post(lambda job: job.isfile("nvt_rdf.png"))
@Project.post(lambda job: job.isfile("nvt_rdf.txt"))
def rdf_nvt_analysis(job):
    """Run RDF analysis."""
    from reproducibility_project.src.analysis.rdf import gsd_rdf

    # RDF
    gsd_rdf(job, ensemble="nvt")


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(
    lambda job: job.doc.get("npt/sampling_results", {}).get("pressure")
)
@flow.with_job
def determine_npt_pressure_sampling(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="npt", prop="pressure", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(
    lambda job: job.doc.get("npt/sampling_results", {}).get("density")
)
@flow.with_job
def determine_npt_density_sampling(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="npt", prop="density", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(
    lambda job: job.doc.get("npt/sampling_results", {}).get("potential_energy")
)
@flow.with_job
def determine_npt_potential_energy_sampling(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="npt", prop="potential_energy", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(
    lambda job: job.doc.get("npt/sampling_results", {}).get("kinetic_energy")
)
@flow.with_job
def determine_npt_kinetic_energy_sampling(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="npt", prop="kinetic_energy", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(
    lambda job: job.doc.get("npt/sampling_results", {}).get("temperature")
)
@flow.with_job
def determine_npt_temperature_sampling(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="npt", prop="temperature", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(
    lambda job: job.doc.get("nvt/sampling_results", {}).get("potential_energy")
)
@flow.with_job
def determine_nvt_potential_energy_sampling(job):
    """Write out sampling results for NVT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="nvt", prop="potential_energy", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(
    lambda job: job.doc.get("nvt/sampling_results", {}).get("kinetic_energy")
)
@flow.with_job
def determine_nvt_kinetic_energy_sampling(job):
    """Write out sampling results for NVT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="nvt", prop="kinetic_energy", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(lambda job: job.doc.get("nvt/sampling_results", {}).get("volume"))
@flow.with_job
def determine_nvt_volume_sampling(job):
    """Write out sampling results for NVT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="nvt", prop="volume", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(
    lambda job: job.doc.get("nvt/sampling_results", {}).get("temperature")
)
@flow.with_job
def determine_nvt_temperature_sampling(job):
    """Write out sampling results for NVT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(
        job=job, ensemble="nvt", prop="temperature", filename=None
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.pre(
    lambda job: all(
        [
            _is_prop_subsampled(job, ensemble="npt", prop=my_prop)
            for my_prop in [
                "potential_energy",
                "kinetic_energy",
                "temperature",
                "pressure",
                "density",
            ]
        ]
    )
)
@Project.post(lambda job: job.doc.get("npt/max_t0"))
def npt_write_largest_t0(job):
    """Write out maximium t0 for all subsampled values in npt production."""
    npt_props = [
        "potential_energy",
        "kinetic_energy",
        "temperature",
        "pressure",
        "density",
    ]
    job.doc["npt/max_t0"] = _get_largest_t0(
        job, ensemble="npt", props=npt_props
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.pre(
    lambda job: all(
        [
            _is_prop_subsampled(job, ensemble="nvt", prop=my_prop)
            for my_prop in [
                "potential_energy",
                "kinetic_energy",
                "temperature",
                "volume",
            ]
        ]
    )
)
@Project.post(lambda job: job.doc.get("nvt/max_t0"))
def nvt_write_largest_t0(job):
    """Write out maximium t0 for all subsampled values in nvt production."""
    nvt_props = ["potential_energy", "kinetic_energy", "temperature", "volume"]
    job.doc["nvt/max_t0"] = _get_largest_t0(
        job, ensemble="nvt", props=nvt_props
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.pre(lambda job: job.doc.get("npt/max_t0"))
@Project.pre(
    lambda job: all(
        [
            _is_prop_subsampled(job, ensemble="npt", prop=my_prop)
            for my_prop in [
                "potential_energy",
                "kinetic_energy",
                "temperature",
                "density",
                "pressure",
            ]
        ]
    )
)
@Project.post(
    lambda job: all(
        [
            job.isfile(fname)
            for fname in [
                "npt_potential_energy.h5",
                "npt_kinetic_energy.h5",
                "npt_temperature.h5",
                "npt_pressure.h5",
                "npt_density.h5",
            ]
        ]
    )
)
@flow.with_job
def npt_write_subsampled_max_t0(job):
    """Write subsampled properties to job.stores."""
    from reproducibility_project.src.analysis.sampler import (
        get_decorr_samples_using_max_t0,
    )

    ensemble = "npt"
    props = [
        "potential_energy",
        "kinetic_energy",
        "temperature",
        "pressure",
        "density",
    ]
    for prop in props:
        with job.stores[f"{ensemble}_{prop}"] as data:
            data["property"] = get_decorr_samples_using_max_t0(
                job,
                ensemble=ensemble,
                property_filename="log-npt.txt",
                property=prop,
            )


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.pre(lambda job: job.doc.get("nvt/max_t0"))
@Project.pre(
    lambda job: all(
        [
            _is_prop_subsampled(job, ensemble="nvt", prop=my_prop)
            for my_prop in [
                "potential_energy",
                "kinetic_energy",
                "temperature",
                "volume",
            ]
        ]
    )
)
@Project.post(
    lambda job: all(
        [
            job.isfile(fname)
            for fname in [
                "nvt_potential_energy.h5",
                "nvt_kinetic_energy.h5",
                "nvt_temperature.h5",
                "nvt_volume.h5",
            ]
        ]
    )
)
@flow.with_job
def nvt_write_subsampled_max_t0(job):
    """Write subsampled properties to job.stores."""
    from reproducibility_project.src.analysis.sampler import (
        get_decorr_samples_using_max_t0,
    )

    ensemble = "nvt"
    props = ["potential_energy", "kinetic_energy", "temperature", "density"]
    for prop in props:
        with job.stores[f"{ensemble}_{prop}"] as data:
            data["property"] = get_decorr_samples_using_max_t0(
                job,
                ensemble=ensemble,
                property_filename="log-nvt.txt",
                property=prop,
            )


@Project.operation
@Project.pre(lambda job: job.isfile("npt_potential_energy.h5"))
@Project.post(lambda job: job.doc.get("npt_potential_energy_avg"))
@Project.post(lambda job: job.doc.get("npt_potential_energy_std"))
@flow.with_job
def npt_calc_potential_energy_statistics(job):
    """Calc statistics on subsampled npt property."""
    _calc_statistics(job, ensemble="npt", prop="potential_energy")


@Project.operation
@Project.pre(lambda job: job.isfile("npt_kinetic_energy.h5"))
@Project.post(lambda job: job.doc.get("npt_kinetic_energy_avg"))
@Project.post(lambda job: job.doc.get("npt_kinetic_energy_std"))
@flow.with_job
def npt_calc_kinetic_energy_statistics(job):
    """Calc statistics on subsampled npt property."""
    _calc_statistics(job, ensemble="npt", prop="kinetic_energy")


@Project.operation
@Project.pre(lambda job: job.isfile("npt_temperature.h5"))
@Project.post(lambda job: job.doc.get("npt_temperature_avg"))
@Project.post(lambda job: job.doc.get("npt_temperature_std"))
@flow.with_job
def npt_calc_temperature_statistics(job):
    """Calc statistics on subsampled npt property."""
    _calc_statistics(job, ensemble="npt", prop="temperature")


@Project.operation
@Project.pre(lambda job: job.isfile("npt_pressure.h5"))
@Project.post(lambda job: job.doc.get("npt_pressure_avg"))
@Project.post(lambda job: job.doc.get("npt_pressure_std"))
@flow.with_job
def npt_calc_pressure_statistics(job):
    """Calc statistics on subsampled npt property."""
    _calc_statistics(job, ensemble="npt", prop="pressure")


@Project.operation
@Project.pre(lambda job: job.isfile("npt_density.h5"))
@Project.post(lambda job: job.doc.get("npt_density_avg"))
@Project.post(lambda job: job.doc.get("npt_density_std"))
@flow.with_job
def npt_calc_density_statistics(job):
    """Calc statistics on subsampled npt property."""
    _calc_statistics(job, ensemble="npt", prop="density")


@Project.operation
@Project.pre(lambda job: job.isfile("nvt_potential_energy.h5"))
@Project.post(lambda job: job.doc.get("nvt_potential_energy_avg"))
@Project.post(lambda job: job.doc.get("nvt_potential_energy_std"))
@flow.with_job
def nvt_calc_potential_energy_statistics(job):
    """Calc statistics on subsampled nvt property."""
    _calc_statistics(job, ensemble="nvt", prop="potential_energy")


@Project.operation
@Project.pre(lambda job: job.isfile("nvt_kinetic_energy.h5"))
@Project.post(lambda job: job.doc.get("nvt_kinetic_energy_avg"))
@Project.post(lambda job: job.doc.get("nvt_kinetic_energy_std"))
@flow.with_job
def nvt_calc_kinetic_energy_statistics(job):
    """Calc statistics on subsampled nvt property."""
    _calc_statistics(job, ensemble="nvt", prop="kinetic_energy")


@Project.operation
@Project.pre(lambda job: job.isfile("nvt_temperature.h5"))
@Project.post(lambda job: job.doc.get("nvt_temperature_avg"))
@Project.post(lambda job: job.doc.get("nvt_temperature_std"))
@flow.with_job
def nvt_calc_temperature_statistics(job):
    """Calc statistics on subsampled nvt property."""
    _calc_statistics(job, ensemble="nvt", prop="temperature")


@Project.operation
@Project.pre(lambda job: job.isfile("nvt_density.h5"))
@Project.post(lambda job: job.doc.get("nvt_density_avg"))
@Project.post(lambda job: job.doc.get("nvt_density_std"))
@flow.with_job
def nvt_calc_density_statistics(job):
    """Calc statistics on subsampled nvt property."""
    _calc_statistics(job, ensemble="nvt", prop="density")


@aggregator.groupby(
    key=[
        "ensemble",
        "engine",
        "molecule",
        "temperature",
        "pressure",
        "cutoff_style",
        "long_range_correction",
    ],
    sort_by="engine",
)
@Project.operation
@Project.pre(
    lambda *jobs: all([job.doc.get("nvt_density_avg") for job in jobs])
)
@Project.pre(
    lambda *jobs: all([job.doc.get("nvt_density_std") for job in jobs])
)
def nvt_calc_aggregate_density_statistics(*jobs):
    """Store aggregate statistics for the nvt density calculations."""
    # should be grouped enough such that only the 16 replicates are the groupings
    agg_project = signac.Project(pathlib.Path("./aggregate_summary/signac.rc"))
    assert (
        len(jobs) == 16
    ), "Not all 16 replicates have their property averaged."
    num_replicas = len(jobs)
    avg_vals = [job.doc.get("nvt_density_avg") for job in jobs]
    avg = np.mean(avg_vals)
    std = np.std(avg_vals)
    sem = std / np.sqrt(num_replicas)
    a_job = jobs[0]
    job_sp = a_job.sp
    _ = job_sp.pop("replica")
    agg_job = agg_project.open_job(statepoint=job_sp)
    agg_job.doc["density-avg"] = avg
    agg_job.doc["density-std"] = std
    agg_job.doc["density-sem"] = sem


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@flow.with_job
def plot_npt_prod_data_with_t0(job):
    """Generate plots for production data with t0 as a vertical line."""
    import pandas as pd

    from reproducibility_project.src.analysis.equilibration import (
        plot_job_property_with_t0,
    )

    ensemble = "npt"

    # plot t0
    df = pd.read_csv(job.fn("log-npt.txt"), delim_whitespace=True, header=0)
    for prop in df.columns:
        data_plt_kwarg = {"label": prop}
        fname = str(prop) + "-" + ensemble + ".png"
        plot_job_property_with_t0(
            job,
            filename=fname,
            property_name=prop,
            log_filename="log-npt.txt",
            title=prop.upper(),
            overwrite=True,
            threshold_fraction=0.0,
            threshold_neff=1,
            vline_scale=1.1,
            data_plt_kwargs=data_plt_kwarg,
        )


@Project.operation
# @Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@flow.with_job
def plot_nvt_prod_data_with_t0(job):
    """Generate plots for production data with t0 as a vertical line."""
    import pandas as pd

    from reproducibility_project.src.analysis.equilibration import (
        plot_job_property_with_t0,
    )

    ensemble = "nvt"

    # plot t0
    df = pd.read_csv(job.fn("log-nvt.txt"), delim_whitespace=True, header=0)
    for prop in df.columns:
        data_plt_kwarg = {"label": prop}
        fname = str(prop) + "-" + ensemble + ".png"
        plot_job_property_with_t0(
            job,
            filename=fname,
            property_name=prop,
            log_filename="log-nvt.txt",
            title=prop.upper(),
            overwrite=True,
            threshold_fraction=0.0,
            threshold_neff=1,
            vline_scale=1.1,
            data_plt_kwargs=data_plt_kwarg,
        )


if __name__ == "__main__":
    pr = Project()
    pr.main()

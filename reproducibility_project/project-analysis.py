"""Setup for signac, signac-flow, signac-dashboard for this study."""
from typing import List

import flow
import numpy as np
import signac


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


def _determine_sampling_information(
    job: signac.contrib.project.Job, ensemble: str, prop: str
) -> None:
    """Write out sampling results for production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    sample_job(
        job,
        ensemble=ensemble,
        filename=f"log-{ensemble}.txt",
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

    _determine_sampling_information(job=job, ensemble="npt", prop="pressure")


@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(
    lambda job: job.doc.get("npt/sampling_results", {}).get("density")
)
@flow.with_job
def determine_npt_density_sampling(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(job=job, ensemble="npt", prop="density")


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
        job=job, ensemble="npt", prop="potential_energy"
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
        job=job, ensemble="npt", prop="kinetic_energy"
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

    _determine_sampling_information(job=job, ensemble="npt", prop="temperature")


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
        job=job, ensemble="nvt", prop="potential_energy"
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
        job=job, ensemble="nvt", prop="kinetic_energy"
    )


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(lambda job: job.doc.get("nvt/sampling_results", {}).get("volume"))
@flow.with_job
def determine_nvt_volume_sampling(job):
    """Write out sampling results for NVT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(job=job, ensemble="nvt", prop="volume")


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(
    lambda job: job.doc.get("nvt/sampling_results", {}).get("temperature")
)
@flow.with_job
def determine_nvt_temperature_sampling(job):
    """Write out sampling results for NVT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(job=job, ensemble="nvt", prop="temperature")


@Project.operation
# @Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
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
# @Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
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
# @Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
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
    props = ["potential_energy", "kinetic_energy", "temperature", "volume"]
    for prop in props:
        with job.stores[f"{ensemble}_{prop}"] as data:
            data["property"] = get_decorr_samples_using_max_t0(
                job,
                ensemble=ensemble,
                property_filename="log-nvt.txt",
                property=prop,
            )


@Project.operation
# @Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
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

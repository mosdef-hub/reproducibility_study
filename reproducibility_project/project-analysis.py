"""Setup for signac, signac-flow, signac-dashboard for this study."""
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


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
def rdf_npt_analysis(job):
    """Run RDF analysis."""
    from reproducibility_project.src.analysis.rdf import gsd_rdf

    # RDF
    gsd_rdf(job)


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
def rdf_nvt_analysis(job):
    """Run RDF analysis."""
    from reproducibility_project.src.analysis.rdf import gsd_rdf

    # RDF
    gsd_rdf(job)


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
@Project.post(lambda job: job.doc.get("npt/sampling_results", {}).get("volume"))
@flow.with_job
def determine_npt_volume_sampling(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    _determine_sampling_information(job=job, ensemble="npt", prop="volume")


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
@Project.pre(lambda job: not job.isfile("pressure.h5"))
@Project.pre(
    lambda job: job.doc.get("npt/sampling_results", {}).get("pressure")
)
@Project.post(lambda job: job.isfile("pressure.h5"))
def write_npt_pressure_subsamples(job):
    """Write pressure subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_pressure as press:
        press["pressure"] = get_subsampled_values(
            job,
            property="pressure",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("volume.h5"))
@Project.pre(lambda job: job.doc.get("npt/sampling_results", {}).get("volume"))
@Project.post(lambda job: job.isfile("volume.h5"))
def write_npt_volume_subsamples(job):
    """Write volume subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_volume as press:
        press["volume"] = get_subsampled_values(
            job,
            property="volume",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("temperature.h5"))
@Project.pre(
    lambda job: job.doc.get("npt/sampling_results", {}).get("temperature")
)
@Project.post(lambda job: job.isfile("temperature.h5"))
def write_npt_temperature_subsamples(job):
    """Write temperature subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_temperature as press:
        press["temperature"] = get_subsampled_values(
            job,
            property="temperature",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("kinetic_energy.h5"))
@Project.pre(
    lambda job: job.doc.get("npt/sampling_results", {}).get("kinetic_energy")
)
@Project.post(lambda job: job.isfile("kinetic_energy.h5"))
def write_npt_kinetic_energy_subsamples(job):
    """Write kinetic_energy subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_kinetic_energy as press:
        press["kinetic_energy"] = get_subsampled_values(
            job,
            property="kinetic_energy",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
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
@Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@flow.with_job
def plot_nvt_prod_data_with_t0(job):
    """Generate plots for production data with t0 as a vertical line."""
    import pandas as pd

    from reproducibility_project.src.analysis.equilibration import (
        plot_job_property_with_t0,
    )

    # plot t0
    df = pd.read_csv(job.fn("log-nvt.txt"), delim_whitespace=True, header=0)
    for prop in df.columns:
        data_plt_kwarg = {"label": prop}
        fname = str(prop) + "-nvt" + ".png"
        plot_job_property_with_t0(
            job,
            filename=fname,
            log_filename="log-nvt.txt",
            property_name=prop,
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

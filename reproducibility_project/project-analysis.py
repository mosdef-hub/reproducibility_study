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


def _is_prop_subsampled(
    job: signac.contrib.Project.Job, ensemble: str, prop: str
):
    """Check if the property has been subsampled."""
    return job.isfile(f"{ensemble}_{prop}.h5")


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
def npt_volume_subsampled(job):
    """Check if volume has been subsampled."""
    return _is_prop_subsampled(job, ensemble="npt", prop="volume")


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
@Project.pre(lambda job: not job.isfile("npt_pressure.h5"))
@Project.pre(
    lambda job: job.doc.get("npt/sampling_results", {}).get("pressure")
)
@Project.post(lambda job: job.isfile("npt_pressure.h5"))
def write_npt_pressure_subsamples(job):
    """Write pressure subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_pressure as fp:
        fp["pressure"] = get_subsampled_values(
            job,
            property="pressure",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("npt_volume.h5"))
@Project.pre(lambda job: job.doc.get("npt/sampling_results", {}).get("volume"))
@Project.post(lambda job: job.isfile("npt_volume.h5"))
def write_npt_volume_subsamples(job):
    """Write volume subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_volume as fp:
        fp["volume"] = get_subsampled_values(
            job,
            property="volume",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("npt_temperature.h5"))
@Project.pre(
    lambda job: job.doc.get("npt/sampling_results", {}).get("temperature")
)
@Project.post(lambda job: job.isfile("npt_temperature.h5"))
def write_npt_temperature_subsamples(job):
    """Write temperature subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_temperature as fp:
        fp["temperature"] = get_subsampled_values(
            job,
            property="temperature",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("npt_kinetic_energy.h5"))
@Project.pre(
    lambda job: job.doc.get("npt/sampling_results", {}).get("kinetic_energy")
)
@Project.post(lambda job: job.isfile("npt_kinetic_energy.h5"))
def write_npt_kinetic_energy_subsamples(job):
    """Write kinetic_energy subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_kinetic_energy as fp:
        fp["kinetic_energy"] = get_subsampled_values(
            job,
            property="kinetic_energy",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("npt_potential_energy.h5"))
@Project.pre(
    lambda job: job.doc.get("npt/sampling_results", {}).get("potential_energy")
)
@Project.post(lambda job: job.isfile("npt_potential_energy.h5"))
def write_npt_potential_energy_subsamples(job):
    """Write potential_energy subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.npt_potential_energy as fp:
        fp["potential_energy"] = get_subsampled_values(
            job,
            property="potential_energy",
            property_filename="log-npt.txt",
            ensemble="npt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("nvt_volume.h5"))
@Project.pre(lambda job: job.doc.get("nvt/sampling_results", {}).get("volume"))
@Project.post(lambda job: job.isfile("nvt_volume.h5"))
def write_nvt_volume_subsamples(job):
    """Write volume subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.nvt_volume as fp:
        fp["volume"] = get_subsampled_values(
            job,
            property="volume",
            property_filename="log-nvt.txt",
            ensemble="nvt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("nvt_temperature.h5"))
@Project.pre(
    lambda job: job.doc.get("nvt/sampling_results", {}).get("temperature")
)
@Project.post(lambda job: job.isfile("nvt_temperature.h5"))
def write_nvt_temperature_subsamples(job):
    """Write temperature subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.nvt_temperature as fp:
        fp["temperature"] = get_subsampled_values(
            job,
            property="temperature",
            property_filename="log-nvt.txt",
            ensemble="nvt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("nvt_kinetic_energy.h5"))
@Project.pre(
    lambda job: job.doc.get("nvt/sampling_results", {}).get("kinetic_energy")
)
@Project.post(lambda job: job.isfile("nvt_kinetic_energy.h5"))
def write_nvt_kinetic_energy_subsamples(job):
    """Write kinetic_energy subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.nvt_kinetic_energy as fp:
        fp["kinetic_energy"] = get_subsampled_values(
            job,
            property="kinetic_energy",
            property_filename="log-nvt.txt",
            ensemble="nvt",
        )


@Project.operation
@Project.pre(lambda job: not job.isfile("nvt_potential_energy.h5"))
@Project.pre(
    lambda job: job.doc.get("nvt/sampling_results", {}).get("potential_energy")
)
@Project.post(lambda job: job.isfile("nvt_potential_energy.h5"))
def write_nvt_potential_energy_subsamples(job):
    """Write kinetic_energy subsamples based on the results of sample_job."""
    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
    )

    with job.stores.nvt_potential_energy as fp:
        fp["potential_energy"] = get_subsampled_values(
            job,
            property="potential_energy",
            property_filename="log-nvt.txt",
            ensemble="nvt",
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

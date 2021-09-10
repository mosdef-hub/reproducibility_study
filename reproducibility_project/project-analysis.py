"""Setup for signac, signac-flow, signac-dashboard for this study."""
import flow
import numpy as np


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


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
@Project.post(lambda job: job.doc["npt/sampling_results"]["volume"])
@Project.post(lambda job: job.doc["npt/sampling_results"]["potential_energy"])
@Project.post(lambda job: job.doc["npt/sampling_results"]["temperature"])
@Project.post(lambda job: job.doc["npt/sampling_results"]["kinetic_energy"])
@Project.post(lambda job: job.doc["npt/sampling_results"]["pressure"])
@flow.with_job
def sample_npt_properties(job):
    """Write out sampling results for NPT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    properties = [
        "volume",
        "potential_energy",
        "temperature",
        "kinetic_energy",
        "pressure",
    ]
    for prop in properties:
        sample_job(
            job,
            ensemble="npt",
            filename="log-npt.txt",
            variable=prop,
            threshold_fraction=0.75,
            threshold_neff=100,
        )


@Project.operation
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(lambda job: job.doc["nvt/sampling_results"]["volume"])
@Project.post(lambda job: job.doc["nvt/sampling_results"]["potential_energy"])
@Project.post(lambda job: job.doc["nvt/sampling_results"]["temperature"])
@Project.post(lambda job: job.doc["nvt/sampling_results"]["kinetic_energy"])
@flow.with_job
def sample_nvt_properties(job):
    """Write out sampling results for NVT production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    properties = ["volume", "potential_energy", "temperature", "kinetic_energy"]
    for prop in properties:
        sample_job(
            job,
            ensemble="nvt",
            filename="log-nvt.txt",
            variable=prop,
            threshold_fraction=0.75,
            threshold_neff=100,
        )


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.post(lambda job: job.data.get("nvt/subsamples/volume", None))
@Project.post(lambda job: job.data.get("nvt/subsamples/potential_energy", None))
@Project.post(lambda job: job.data.get("nvt/subsamples/temperature", None))
@Project.post(lambda job: job.data.get("nvt/subsamples/kinetic_energy", None))
@flow.with_job
def write_nvt_properties(job):
    """Write subsampled nvt property data points to the job.data store."""
    from reproducibility_project.src.analysis.sampler import sample_job

    properties = ["potential_energy", "temperature", "volume", "kinetic_energy"]
    for prop in properties:
        sample_job(
            job,
            filename="log-nvt.txt",
            variable=prop,
            threshold_fraction=0.75,
            threshold_neff=100,
        )


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(lambda job: job.data.get("npt/subsamples/volume", None))
@Project.post(lambda job: job.data.get("npt/subsamples/potential_energy", None))
@Project.post(lambda job: job.data.get("npt/subsamples/temperature", None))
@Project.post(lambda job: job.data.get("npt/subsamples/kinetic_energy", None))
@Project.post(lambda job: job.data.get("npt/subsamples/pressure", None))
@flow.with_job
def write_npt_properties(job):
    """Write subsampled npt property data points to the job.data store."""
    from reproducibility_project.src.analysis.sampler import sample_job

    properties = [
        "potential_energy",
        "temperature",
        "volume",
        "kinetic_energy",
        "pressure",
    ]
    for prop in properties:
        sample_job(
            job,
            filename="log-npt.txt",
            variable=prop,
            threshold_fraction=0.75,
            threshold_neff=100,
        )


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.pre(lambda job: job.data.get("nvt/subsamples/volume", None))
@Project.pre(lambda job: job.data.get("nvt/subsamples/potential_energy", None))
@Project.pre(lambda job: job.data.get("nvt/subsamples/temperature", None))
@Project.pre(lambda job: job.data.get("nvt/subsamples/kinetic_energy", None))
@Project.post(lambda job: job.doc.sp.get("volume-nvt-avg", None))
@Project.post(lambda job: job.doc.sp.get("volume-nvt-std", None))
@Project.post(lambda job: job.doc.sp.get("potential_energy-nvt-avg", None))
@Project.post(lambda job: job.doc.sp.get("potential_energy-nvt-std", None))
@Project.post(lambda job: job.doc.sp.get("temperature-nvt-avg", None))
@Project.post(lambda job: job.doc.sp.get("temperature-nvt-std", None))
@Project.post(lambda job: job.doc.sp.get("kinetic_energy-nvt-avg", None))
@Project.post(lambda job: job.doc.sp.get("kinetic_energy-nvt-std", None))
@flow.with_job
def average_nvt_sampled_data(job):
    """Write avg and std of the sampled job.data values."""
    properties = ["potential_energy", "temperature", "volume", "kinetic_energy"]
    for prop in properties:
        with job.data:
            job.doc[f"{prop}-nvt-avg"] = np.mean(
                job.data[f"subsamples/{prop}-nvt"]
            )
            job.doc[f"{prop}-nvt-std"] = np.std(
                job.data[f"subsamples/{prop}-nvt"]
            )


@Project.operation
@Project.pre(lambda job: job.isfile("trajectory-npt.gsd"))
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.pre(lambda job: job.data.get("npt/subsamples/volume", None))
@Project.pre(lambda job: job.data.get("npt/subsamples/potential_energy", None))
@Project.pre(lambda job: job.data.get("npt/subsamples/temperature", None))
@Project.pre(lambda job: job.data.get("npt/subsamples/kinetic_energy", None))
@Project.pre(lambda job: job.data.get("npt/subsamples/pressure", None))
@Project.post(lambda job: job.doc.sp.get("volume-npt-avg", None))
@Project.post(lambda job: job.doc.sp.get("volume-npt-std", None))
@Project.post(lambda job: job.doc.sp.get("potential_energy-npt-avg", None))
@Project.post(lambda job: job.doc.sp.get("potential_energy-npt-std", None))
@Project.post(lambda job: job.doc.sp.get("temperature-npt-avg", None))
@Project.post(lambda job: job.doc.sp.get("temperature-npt-std", None))
@Project.post(lambda job: job.doc.sp.get("kinetic_energy-npt-avg", None))
@Project.post(lambda job: job.doc.sp.get("kinetic_energy-npt-std", None))
@Project.post(lambda job: job.doc.sp.get("pressure-npt-avg", None))
@Project.post(lambda job: job.doc.sp.get("pressure-npt-std", None))
def average_npt_sampled_data(job):
    """Write avg and std of the sampled job.data values."""
    properties = [
        "potential_energy",
        "temperature",
        "volume",
        "kinetic_energy",
        "pressure",
    ]
    for prop in properties:
        with job.data:
            job.doc[f"{prop}-npt-avg"] = np.mean(
                job.data[f"subsamples/{prop}-npt"]
            )
            job.doc[f"{prop}-npt-std"] = np.std(
                job.data[f"subsamples/{prop}-npt"]
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

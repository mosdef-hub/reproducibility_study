"""Setup for signac, signac-flow, signac-dashboard for this study."""

import pathlib
from copy import copy
from typing import List

import flow
import numpy as np
import signac
from flow import aggregator
from flow.environment import DefaultPBSEnvironment


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


class RahmanAnalysis(DefaultPBSEnvironment):
    """Subclass of DefaultPBSEnvironment for VU's Rahman cluster."""

    hostname_pattern = "master.cl.vanderbilt.edu"
    template = "rahman_analysis.sh"

    @classmethod
    def add_args(cls, parser):
        """Add command line arguments to the submit call."""
        parser.add_argument(
            "--walltime",
            type=float,
            default=96,
            help="Walltime for this submission",
        )


mc_engines = ("gomc", "mcccs", "cassandra")
md_engines = ("gromacs", "hoomd", "lammps-VU")  # , "lammps-UD")
md_npt_props = [
    "potential_energy",
    "kinetic_energy",
    "temperature",
    "pressure",
    "density",
]
mc_npt_props = [
    "potential_energy",
    "temperature",
    "pressure",
    "density",
]
# md_nvt_props = [
#     "potential_energy",
#     "kinetic_energy",
#     "temperature",
#     # "density",
# ]
# mc_nvt_props = [
#     "potential_energy",
#     "temperature",
#     # "density",
# ]


# make a FlowGroup for all analysis operations, usually faster than the plotting methods
analysis = Project.make_group(name="analysis")


def _aggregate_statistics(val, n):
    avg = np.mean(val)
    std = np.std(val)
    sem = std / np.sqrt(n)
    pass


def _determine_sampling_information(
    job: signac.contrib.project.Job,
    ensemble: str,
    prop: str,
    strict: bool,
    filename: str = None,
) -> None:
    """Write out sampling results for production properties."""
    from reproducibility_project.src.analysis.sampler import sample_job

    if filename is None:
        filename = f"log-{ensemble}.txt"
    # Monte Carlo groups used pyMBAR to convert their log-$ENSEMBLE.txt files to production data already
    # The sampling rates and lack of time-dependency makes pyMBAR more difficult to use for these prodution data
    # see sample_job for information about this
    if job.sp.engine in mc_engines:
        sample_job(
            job,
            ensemble=ensemble,
            filename=filename,
            variable=prop,
            threshold_fraction=0.75,
            threshold_neff=100,
            strict=strict,
            monte_carlo_override=True,
        )
    else:
        sample_job(
            job,
            ensemble=ensemble,
            filename=filename,
            variable=prop,
            threshold_fraction=0.75,
            threshold_neff=100,
            strict=strict,
            monte_carlo_override=False,
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
    if job.sp.engine in mc_engines:
        props = mc_npt_props
    else:
        props = md_npt_props
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


# @Project.label
# def nvt_ke_subsampled(job):
#     """Check if kinetic energy has been subsampled."""
#     return _is_prop_subsampled(job, ensemble="nvt", prop="kinetic_energy")


@Project.label
def npt_pe_subsampled(job):
    """Check if potential energy has been subsampled."""
    return _is_prop_subsampled(job, ensemble="npt", prop="potential_energy")

    # @Project.label
    # def nvt_pe_subsampled(job):
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

    # @Project.label
    # def nvt_temp_subsampled(job):
    """Check if temperature has been subsampled."""
    return _is_prop_subsampled(job, ensemble="nvt", prop="temperature")

    # @Project.label
    # def nvt_volume_subsampled(job):
    """Check if volume has been subsampled."""
    return _is_prop_subsampled(job, ensemble="nvt", prop="volume")


@Project.label
def npt_prod_finished(job):
    """Generate label if npt production is complete."""
    return job.isfile("trajectory-npt.gsd") and job.isfile("log-npt.txt")

    # @Project.label
    # def nvt_prod_finished(job):
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

    # @Project.label
    # def nvt_rdf_calculated(job):
    """Generate label if nvt rdf calc is complete."""
    return (
        job.isfile("trajectory-nvt.gsd")
        and job.isfile("nvt_rdf.png")
        and job.isfile("nvt_rdf.txt")
    )


def create_rdf_analysis(ensemble: str):
    """Dynamically create the rdf analysis methods."""

    @Project.operation(f"rdf-{ensemble}")
    @Project.pre(lambda job: job.isfile(f"trajectory-{ensemble}.gsd"))
    @Project.post(lambda job: job.isfile(f"{ensemble}_rdf.png"))
    @Project.post(lambda job: job.isfile(f"{ensemble}_rdf.txt"))
    def rdf_analysis(job):
        """Run RDF analysis."""
        from reproducibility_project.src.analysis.rdf import gsd_rdf

        # RDF
        gsd_rdf(job, filename=f"trajectory-{ensemble}.gsd", ensemble=ensemble)


for ensemble in ["npt", "nvt"]:
    create_rdf_analysis(ensemble=ensemble)


def create_property_sampling(
    ensemble: str,
    prop: str,
    simulation_type: str,
    engine_list: List[str],
):
    """Dynamically create the sampling steps for the simulation data."""

    @analysis
    @Project.operation(
        f"{simulation_type}_determine_{ensemble}_{prop}_sampling"
    )
    @Project.pre(lambda job: job.sp.engine in engine_list)
    @Project.pre(lambda job: job.isfile(f"log-{ensemble}.txt"))
    @Project.post(
        lambda job: job.doc.get(f"{ensemble}/sampling_results", {}).get(
            f"{prop}"
        )
        is not None
    )
    @flow.with_job
    def sample_property(job):
        """Write out sampling results for production properties."""
        from reproducibility_project.src.analysis.sampler import sample_job

        _determine_sampling_information(
            job=job,
            ensemble=ensemble,
            prop=prop,
            filename=None,
            strict=False,
        )


# generate md sampling methods
# for ensemble, props in zip(["npt", "nvt"], [md_npt_props, md_nvt_props]):
for ensemble, props in zip(["npt"], [md_npt_props]):
    for prop in props:
        create_property_sampling(
            ensemble=ensemble,
            prop=prop,
            simulation_type="md",
            engine_list=md_engines,
        )

# generate mc sampling methods
# for ensemble, props in zip(["npt", "nvt"], [mc_npt_props, mc_nvt_props]):
for ensemble, props in zip(["npt"], [mc_npt_props]):
    for prop in props:
        create_property_sampling(
            ensemble=ensemble,
            prop=prop,
            simulation_type="mc",
            engine_list=mc_engines,
        )


def create_largest_t0_operations(
    ensemble: str,
    prop_list: List[str],
    simulation_type: str,
    engine_list: List[str],
):
    """Dynamically create the operations to determine the largest t0 for sampling."""

    @analysis
    @Project.operation(f"{ensemble}_{simulation_type}_write_largest_t0")
    @Project.pre(lambda job: job.isfile(f"log-{ensemble}.txt"))
    @Project.pre(lambda job: job.sp.engine in engine_list)
    @Project.pre(
        lambda job: all(
            [
                _is_prop_subsampled(job, ensemble=ensemble, prop=my_prop)
                for my_prop in prop_list
            ]
        )
    )
    @Project.post(lambda job: job.doc.get(f"{ensemble}/max_t0") is not None)
    def write_largest_t0(job):
        """Write out maximium t0 for all subsampled values in npt production."""
        job.doc[f"{ensemble}/max_t0"] = _get_largest_t0(
            job, ensemble=ensemble, props=prop_list
        )


# md operations
# for ensemble, prop_list in zip(["npt", "nvt"], [md_npt_props, md_nvt_props]):
for ensemble, prop_list in zip(["npt"], [md_npt_props]):
    create_largest_t0_operations(
        ensemble=ensemble,
        prop_list=prop_list,
        simulation_type="md",
        engine_list=md_engines,
    )

# mc operations
# for ensemble, prop_list in zip(["npt", "nvt"], [mc_npt_props, mc_nvt_props]):
for ensemble, prop_list in zip(["npt"], [mc_npt_props]):
    create_largest_t0_operations(
        ensemble=ensemble,
        prop_list=prop_list,
        simulation_type="mc",
        engine_list=mc_engines,
    )


def create_write_subsampled_max_t0(
    ensemble: str,
    prop_list: List[str],
    simulation_type: str,
    engine_list: List[str],
    is_monte_carlo: bool,
):
    """Dynamically create the operations to write the subsampled data based on max t0."""

    @analysis
    @Project.operation(
        f"{ensemble}_{simulation_type}_write_subsampled_data_max_t0"
    )
    @Project.pre(lambda job: job.isfile(f"log-{ensemble}.txt"))
    @Project.pre(lambda job: job.doc.get(f"{ensemble}/max_t0") is not None)
    @Project.pre(lambda job: job.sp.engine in engine_list)
    @Project.post(
        lambda job: all(
            [job.isfile(f"{ensemble}_{prop}.h5") for prop in prop_list]
        )
    )
    @flow.with_job
    def write_subsampled_data_max_t0(job):
        """Write subsampled properties to job.stores."""
        from reproducibility_project.src.analysis.sampler import (
            get_decorr_samples_using_max_t0,
        )

        for prop in prop_list:
            with job.stores[f"{ensemble}_{prop}"] as data:
                data["property"] = get_decorr_samples_using_max_t0(
                    job,
                    ensemble=ensemble,
                    property_filename=f"log-{ensemble}.txt",
                    prop=prop,
                    is_monte_carlo=is_monte_carlo,
                )


# md operations
# for ensemble, prop_list in zip(["npt", "nvt"], [md_npt_props, md_nvt_props]):
for ensemble, prop_list in zip(["npt"], [md_npt_props]):
    create_write_subsampled_max_t0(
        ensemble=ensemble,
        prop_list=prop_list,
        simulation_type="md",
        engine_list=md_engines,
        is_monte_carlo=False,
    )

# mc operations
# for ensemble, prop_list in zip(["npt", "nvt"], [mc_npt_props, mc_nvt_props]):
for ensemble, prop_list in zip(["npt"], [mc_npt_props]):
    create_write_subsampled_max_t0(
        ensemble=ensemble,
        prop_list=prop_list,
        simulation_type="mc",
        engine_list=mc_engines,
        is_monte_carlo=True,
    )


def create_calc_prop_statistics(
    ensemble: str, prop: str, simulation_type: str, engine_list: List[str]
):
    """Dynamically create the functions to calculate property statistics."""

    @analysis
    @Project.operation(
        f"{simulation_type}_{ensemble}_{prop}_calc_prop_statistics"
    )
    @Project.pre(lambda job: job.isfile(f"{ensemble}_{prop}.h5"))
    @Project.pre(lambda job: job.sp.engine in engine_list)
    @Project.post(lambda job: job.doc.get(f"{ensemble}_{prop}_avg") is not None)
    @Project.post(lambda job: job.doc.get(f"{ensemble}_{prop}_std") is not None)
    @flow.with_job
    def calc_prop_stats(job):
        """Calc statistics on subsampled property."""
        _calc_statistics(job=job, ensemble=ensemble, prop=prop)


# generate md property statistics methods
# for ensemble, props in zip(["npt", "nvt"], [md_npt_props, md_nvt_props]):
for ensemble, prop_list in zip(["npt"], [md_npt_props]):
    for prop in props:
        create_calc_prop_statistics(
            ensemble=ensemble,
            prop=prop,
            simulation_type="md",
            engine_list=md_engines,
        )


# generate mc property statistics methods
# for ensemble, props in zip(["npt", "nvt"], [mc_npt_props, mc_nvt_props]):
for ensemble, prop_list in zip(["npt"], [mc_npt_props]):
    for prop in props:
        create_calc_prop_statistics(
            ensemble=ensemble,
            prop=prop,
            simulation_type="mc",
            engine_list=mc_engines,
        )


def create_calc_aggregate_statistics(
    ensemble: str,
    prop: str,
    simulation_type: str,
):
    """Dynamically create the operations to calculate the aggregate property statistics for each engine."""

    @analysis
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
    @Project.operation(
        f"{ensemble}_{prop}_{simulation_type}_calc_aggregate_stats"
    )
    @Project.pre(
        lambda *jobs: all(
            [job.doc.get(f"{ensemble}_{prop}_avg") is not None for job in jobs]
        )
    )
    @Project.pre(
        lambda *jobs: all(
            [job.doc.get(f"{ensemble}_{prop}_std") is not None for job in jobs]
        )
    )
    @Project.post(
        lambda *jobs: all(
            [job.isfile(f"{ensemble}_{prop}_calculated.txt") for job in jobs]
        )
    )
    def calc_aggregate_prop_statistics(*jobs):
        """Store aggregate statistics for the ensemble, property calculations."""
        # should be grouped enough such that only the 16 replicates are the groupings
        agg_project = signac.get_project(
            root=pathlib.Path("./aggregate_summary/").absolute()
        )
        assert (
            len(jobs) == 16
        ), "Not all 16 replicates have their property averaged."
        num_replicas = len(jobs)
        avg_vals = [job.doc.get(f"{ensemble}_{prop}_avg") for job in jobs]
        avg = np.mean(avg_vals)
        std = np.std(avg_vals)
        sem = std / np.sqrt(num_replicas)
        a_job = jobs[0]
        job_sp = copy(a_job.sp._to_base())  # convert from syncedDict->dict

        # Discard unnesscary properties
        discarded_replica = job_sp.pop("replica")
        discarded_mass = job_sp.pop("mass")

        agg_job = agg_project.open_job(statepoint=job_sp)
        agg_job.doc[f"{prop}-avg"] = avg
        agg_job.doc[f"{prop}-std"] = std
        agg_job.doc[f"{prop}-sem"] = sem
        for job in jobs:
            my_job_path = pathlib.Path(
                job.fn(f"{ensemble}_{prop}_calculated.txt")
            ).touch()


# generate md aggreagate values
# for ensemble, prop_list in zip(["npt", "nvt"], [md_npt_props, md_nvt_props]):
for ensemble, prop_list in zip(["npt"], [md_npt_props]):
    for prop in prop_list:
        create_calc_aggregate_statistics(
            ensemble=ensemble,
            prop=prop,
            simulation_type="md",
        )

# generate mc aggreagate values
# for ensemble, prop_list in zip(["npt", "nvt"], [mc_npt_props, mc_nvt_props]):
for ensemble, prop_list in zip(["npt"], [mc_npt_props]):
    for prop in prop_list:
        create_calc_aggregate_statistics(
            ensemble=ensemble,
            prop=prop,
            simulation_type="mc",
        )

'''
@Project.operation
@Project.pre(lambda job: job.isfile("log-npt.txt"))
@Project.post(lambda job: job.isfile("density-npt.png"))
@flow.with_job
def plot_npt_prod_data_with_t0(job):
    """Generate plots for production data with t0 as a vertical line."""
    import pandas as pd

    from reproducibility_project.src.analysis.equilibration import (
        plot_job_property_with_t0,
    )

    ensemble = "npt"

    # plot t0
    with open(job.fn("log-npt.txt"), "r") as f:
        line1 = f.readline()
        df = pd.read_csv(
            f, delim_whitespace=True, names=line1.replace("#", "").split()
        )
    """
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
            strict=False,
            vline_scale=1.1,
            data_plt_kwargs=data_plt_kwarg,
        )
    """
    data_plt_kwarg = {"label": "density"}
    fname = "density" + "-" + ensemble + ".png"
    plot_job_property_with_t0(
        job,
        filename=fname,
        property_name="density",
        log_filename="log-npt.txt",
        title="Density",
        overwrite=True,
        threshold_fraction=0.0,
        threshold_neff=1,
        strict=False,
        vline_scale=1.1,
        data_plt_kwargs=data_plt_kwarg,
    )
'''

'''
@Project.operation
# @Project.pre(lambda job: job.isfile("trajectory-nvt.gsd"))
@Project.pre(lambda job: job.isfile("log-nvt.txt"))
@Project.pre(lambda job: False)
@flow.with_job
def plot_nvt_prod_data_with_t0(job):
    """Generate plots for production data with t0 as a vertical line."""
    import pandas as pd

    from reproducibility_project.src.analysis.equilibration import (
        plot_job_property_with_t0,
    )

    ensemble = "nvt"

    # plot t0
    with open(job.fn("log-nvt.txt", 'r') as f:
        line1 = f.readline()
        df = pd.read_csv(f, delim_whitespace=True, names=line1.replace('#', '').split())
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
            strict=False,
            vline_scale=1.1,
            data_plt_kwargs=data_plt_kwarg,
        )
'''

if __name__ == "__main__":
    pr = Project()
    pr.main()

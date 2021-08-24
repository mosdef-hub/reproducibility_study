"""Setup for signac, signac-flow, signac-dashboard for this study."""
import flow


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


@Project.operation
@Project.pre(lambda j: j.isfile("trajectory.gsd"))
def rdf_analysis(job):
    """Run analysis."""
    from reproducibility_project.src.analysis.rdf import gsd_rdf

    # RDF
    gsd_rdf(job)


@Project.operation
@Project.pre(lambda j: j.isfile("trajectory.gsd"))
@Project.pre(lambda j: j.isfile("log.txt"))
def plot_prod_data_with_t0(job):
    """Generate plots for production data with t0 as a vertical line."""
    import pandas as pd

    from reproducibility_project.src.analysis.equilibration import (
        plot_job_property_with_t0,
    )

    # plot t0
    df = pd.read_csv(job.fn("log.txt"), delim_whitespace=True, header=0)
    for prop in df.columns:
        data_plt_kwarg = {"label": prop}
        fname = str(prop) + ".png"
        plot_job_property_with_t0(
            job,
            filename=fname,
            property_name=prop,
            title=prop.upper(),
            overwrite=True,
            threshold=0.0,
            vline_scale=1.5,
            data_plt_kwargs=data_plt_kwarg,
        )


if __name__ == "__main__":
    pr = Project()
    pr.main()

"""Setup for signac, signac-flow, signac-dashboard for this study."""
import flow


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


@Project.operation
@Project.pre(lambda j: j.isfile("trajectory.gsd"))
def run_analysis(job):
    """Run analysis."""
    from reproducibility_project.src.analysis.rdf import gsd_rdf

    # RDF
    gsd_rdf(job)


if __name__ == "__main__":
    pr = Project()
    pr.main()

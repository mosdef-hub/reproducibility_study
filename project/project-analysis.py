"""Setup for signac, signac-flow, signac-dashboard for this study."""
from os import path

import flow


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


@Project.operation
@Project.pre(lambda j: path.isfile(j.fn("trajectory.gsd")))
def run_analysis(job):
    """Run analysis."""
    from src.analysis.rdf import gsd_rdf

    # RDF
    pass


if __name__ == "__main__":
    pr = Project()
    pr.main()

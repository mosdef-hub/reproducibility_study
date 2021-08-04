"""Setup for signac, signac-flow, signac-dashboard for this study."""
# import foyer
import os
import pathlib

import flow
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()

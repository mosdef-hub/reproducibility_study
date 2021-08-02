"""Setup for signac, signac-flow, signac-dashboard for this study."""
import flow
from flow import environments
#import foyer
import os
import pathlib


class Project(flow.FlowProject):
    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / 'data'
        self.ff_fn = self.data_dir / 'forcefield.xml'


if __name__ == '__main__':
    pr = Project()
    pr.main()
    breakpoint()

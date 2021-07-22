"""Setup for signac, signac-flow, signac-dashboard for this study."""
import flow
from flow import FlowProject, directives, cmd
from flow import environments


class Project(FlowProject):
    pass


if __name__ == '__main__':
    pr = Project()
    pr.main()

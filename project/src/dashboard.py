"""Signac project visualization support using signac-dashboard."""
from signac_dashboard import Dashboard
from signac_dashboard.modules import (
    ImageViewer,
    Notes,
    StatepointList,
    TextDisplay,
    VideoViewer,
)


class PlotDashboard(Dashboard):
    """Take signac dashboard and plot statepoints."""

    def job_sorter(self, job):
        """Accumulate all jobs according to a statepoint."""
        return [job.sp.simulation_engine]


if __name__ == "__main__":
    modules = []
    modules.append(StatepointList())
    PlotDashboard(modules=modules).main()

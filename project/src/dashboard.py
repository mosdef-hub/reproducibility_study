"""Dashboard for viewing data associated with jobs in a signac project."""
from signac_dashboard import Dashboard
from signac_dashboard.modules import (
    ImageViewer,
    Notes,
    StatepointList,
    TextDisplay,
    VideoViewer,
)


class PlotDashboard(Dashboard):
    """Subclass of main Dashboard class to provide project-specific methods."""

    def job_sorter(self, job):
        """Sorts jobs based on provided statepoint parameters."""
        return [job.sp.simulation_engine]


if __name__ == "__main__":
    modules = []
    modules.append(StatepointList())
    PlotDashboard(modules=modules).main()

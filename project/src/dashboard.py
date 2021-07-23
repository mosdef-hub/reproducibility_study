from signac_dashboard import Dashboard
from signac_dashboard.modules import StatepointList, ImageViewer, VideoViewer
from signac_dashboard.modules import Notes, TextDisplay


class PlotDashboard(Dashboard):
    def job_sorter(self, job):
        return [job.sp.simulation_engine]


if __name__ == '__main__':
    modules = []
    modules.append(StatepointList())
    PlotDashboard(modules=modules).main()

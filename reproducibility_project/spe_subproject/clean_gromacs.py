"""Reset gromacs job."""
import signac

project = signac.get_project()
for job in project.find_jobs({"engine": "gromacs"}):
    print(job)
    job.reset()

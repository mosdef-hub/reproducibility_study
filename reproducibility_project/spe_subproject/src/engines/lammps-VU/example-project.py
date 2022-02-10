"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib

import flow
import numpy as np
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


# ____________________________________________________________________________
"""Setting progress label"""


@Project.label
@Project.pre(lambda j: j.sp.engine == "MYENGINENAME")
def CreatedEngineInput(job):
    """Check if the .json molecule topology was converted to engine input."""
    return job.isfile("MYENGINEINPUTS")

@Project.label
@Project.pre(lambda j: j.sp.engine == "MYENGINENAME")
def OutputThermoData(job):
    """Check if the engine loaded the input files and wrote out thermo data."""
    return job.isfile("NAMEOFMYTHERMOOUTPUT")

@Project.label
@Project.pre(lambda j: j.sp.engine == "MYENGINENAME")
def FinishedSPECalc(job):
    """Check if the log-spe.txt has been created."""
    return job.isfile("log-spe.txt")


# _____________________________________________________________________
"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.engine == "MYENGINENAME")
@Project.post(CreatedEngineInput)
@flow.with_job
def LoadSystemSnapShot(job):
    """Create initial configurations of the system statepoint."""
    import mbuild as mb

    pr = Project()
    snapshot_directory = pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    molecule = job.sp.molecule
    molecule_filename = molecule + '.json'
    box = mb.load(str(snapshot_directory / molecule_filename))
    # Apply forcefield and write out engine input files
    #__________________________________________________

@Project.operation
@Project.pre(lambda j: j.sp.engine == "MYENGINENAME")
@Project.pre(CreatedEngineInput)
@Project.post(OutputThermoData)
@flow.with_job
def CalculateEnergy(job):
    """Load onto a cluster and output the point energy for the snapshot."""
    #__________________________________________________

@Project.operation
@Project.pre(lambda j: j.sp.engine == "MYENGINENAME")
@Project.pre(OutputThermoData)
@Project.post(FinishedSPECalc)
@flow.with_job
@flow.cmd
def FormatTextFile(job):
    """Take the output from the simulation engine and convert it to log-spe.txt for data comparisons.

    See README.md for spe_subproject for formatting information.
    """
    #__________________________________________________

if __name__ == "__main__":
    pr = Project()
    pr.main()

"""Make MCCCS-MN restart file for calculating pressure."""

import mbuild as mb
from fort77maker_onebox import fort77writer
from system_builder import get_molecule

box_length = 1.23379050108225  # nm
initial_filename = "box1config1a.xyz"

fort77writer(
    [get_molecule("methaneUA")],
    box_length,
    mb.load(initial_filename),
)

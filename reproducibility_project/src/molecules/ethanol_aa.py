"""Create atomistic representation of ethanol."""

import os

import mbuild as mb

from reproducibility_project.src import molecules


class EthanolAA(mb.Compound):
    """Create a single atomistic ethanol compound."""

    def __init__(self):
        super(EthanolAA, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/ethanol_aa.mol2"), label="ETO")


def main():
    """Create a EthanolAA compound and print basic properties."""
    ethanol = EthanolAA()
    print(ethanol)
    print(ethanol.name)
    print(ethanol.labels)
    print(ethanol["ETO"])


if __name__ == "__main__":
    main()

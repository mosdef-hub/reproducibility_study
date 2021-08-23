"""Create united-atom representation of benzene."""

import os

import mbuild as mb

from reproducibility_project.src import molecules


class BenzeneUA(mb.Compound):
    """Create a single benzene compound."""

    def __init__(self):
        super(BenzeneUA, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/benzene_ua.mol2"))


def main():
    """Create a BenzeneUA compound and print basic properties."""
    benzene = BenzeneUA()
    print(benzene)
    print(benzene.name)
    print(benzene.labels)
    print(benzene["BEN"])


if __name__ == "__main__":
    main()

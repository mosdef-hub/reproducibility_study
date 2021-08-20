"""Create united-atom representation of benzene."""

import mbuild as mb


class BenzeneUA(mb.Compound):
    """Create a single benzene compound."""

    def __init__(self):
        super(BenzeneUA, self).__init__()
        self.add(mb.load("benzene_ua.mol2"), label="BENZ")


def main():
    """Create a BenzeneUA compound and print basic properties."""
    benzene = BenzeneUA()
    print(benzene)
    print(benzene.name)
    print(benzene.labels)
    print(benzene["BENZ"])


if __name__ == "__main__":
    main()

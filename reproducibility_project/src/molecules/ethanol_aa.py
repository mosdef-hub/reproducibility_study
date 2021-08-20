"""Create atomistic representation of ethanol."""

import mbuild as mb


class EthanolAA(mb.Compound):
    """Create a single atomistic ethanol compound."""

    def __init__(self):
        super(EthanolAA, self).__init__()
        self.add(mb.load("ethanol_aa.mol2", labels="ETO"))


def main():
    """Create a EthanolAA compound and print basic properties."""
    ethanol = EthanolAA()
    print(ethanol)
    print(ethanol.name)
    print(ethanol.labels)
    print(ethanol["ETO"])


if __name__ == "__main__":
    main()

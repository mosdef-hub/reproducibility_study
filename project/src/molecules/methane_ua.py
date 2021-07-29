"""Create united-atom representation of methane."""

import mbuild as mb


class MethaneUA(mb.Compound):
    """Create a single particle methane compound."""

    def __init__(self):
        super(MethaneUA, self).__init__()
        methane = mb.Compound(
            name="_CH4",
        )
        self.add(methane, label="MET")


def main():
    """Create a MethaneUA compound and print basic properties."""
    methane = MethaneUA()
    print(methane)
    print(methane.name)
    print(methane.labels)
    print(methane["MET"])


if __name__ == "__main__":
    main()

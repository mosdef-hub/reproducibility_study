"""Create united-atom representation of pentane."""

import mbuild as mb
import numpy as np


class PentaneUA(mb.Compound):
    """Create a single particle pentane compound."""

    def __init__(self):
        super(PentaneUA, self).__init__()
        # Calculate the angle between the two ports
        angle = np.deg2rad(114)
        x = 0
        y = 0.077
        z = 0
        vec = [
            x,
            y * np.cos(angle) - z * np.sin(angle),
            y * np.sin(angle) + z * np.cos(angle),
        ]

        # Create the end group compound
        ch3 = mb.Compound()
        ch3.add(mb.Particle(name="_CH3"))
        ch3.add(mb.Port(anchor=ch3[0]), "up")
        ch3["up"].translate([x, y, z])

        # Create the internal monomer
        ch2 = mb.Compound()
        ch2.add(mb.Particle(name="_CH2"))
        ch2.add(mb.Port(anchor=ch2[0]), "up")
        ch2["up"].translate([x, y, z])
        ch2.add(mb.Port(anchor=ch2[0], orientation=vec), "down")
        ch2["down"].translate(vec)

        pentane = mb.recipes.Polymer(
            monomers=[ch2], end_groups=[ch3, mb.clone(ch3)]
        )
        pentane.build(n=3)

        self.add(pentane, label="PNT")


def main():
    """Create a PentaneUA compound and print basic properties."""
    pentane = PentaneUA()
    print(pentane)
    print(pentane.name)
    print(pentane.labels)
    print(pentane["PNT"])


if __name__ == "__main__":
    main()

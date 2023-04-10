"""Takes an mbuild filled_box and creates a restart file for mcccs simulations."""
import errno
import os
import sys
from collections import Counter
from subprocess import call

import mbuild as mb
import numpy as np

from reproducibility_project.src.engine_input.mcccs.utils.fort77helpfun import (
    unique,
    xyzwriter,
)


def fort77writer(
    molecules,
    filled_boxes,
    output_file="config.new",
    NBox=2,
    xyz_file=["initial_structure_box1.xyz", "initial_structure_box2.xyz"],
):
    """Create a fort.77 (restart file for MCCCS-MN code) from a structure.

    Parameters
    ----------
    molecules    : List of mBuild molecules
    filled_boxes : List of mBuild filled boxes
    output_file  : Name of the fort.77 file created
    NBox         : Number of boxes in the simulation, currently only supports two box
    xyz_file     : List of the names of the xyz files in which coordinates of the filled_box are stored
    """
    # preparing name_list (or element list from the names of the filled box)
    name_list = []
    for particle in filled_boxes[0].particles():
        name_list.append(particle.name)
    initial_coord = filled_boxes[0].xyz * 10
    xyzwriter(name_list, initial_coord, xyz_file[0])
    name_list = []
    for particle in filled_boxes[1].particles():
        name_list.append(particle.name)
    initial_coord = filled_boxes[1].xyz * 10
    xyzwriter(name_list, initial_coord, xyz_file[1])
    # high-precision xyz files written

    lengths1 = filled_boxes[0].box.lengths
    lengths2 = filled_boxes[1].box.lengths
    mols = []
    for mol in filled_boxes[0].children:
        mols.append(mol.name)
    num_each_moltype_box1 = list(Counter(mols).values())
    mols = []
    for mol in filled_boxes[1].children:
        mols.append(mol.name)
    num_each_moltype_box2 = list(Counter(mols).values())

    #######
    lortho = [True, True]  # nbox elements
    lconfig_file = [
        False,
        False,
    ]  ## are you providing any config file for your box, config file is written by MCCCS-MN
    CellVec = []
    CellLengths = []

    for i in lortho:
        CellVec.append(np.identity(3))

    for i in range(NBox):
        if lortho[i]:
            CellLengths.append(np.ones((1, 3)))
        else:
            CellLengths.append(np.identity(3))
    ###########Change code after this line
    # you need to input cellvecs (3x3 matrix) if your box is a non ortho box
    # if your box is non ortho, input your celllengths as 3x3 matrix, else celllengths[i] would be 1x3 np matrix
    # in this example cellvec[0] is a 3x3 matrix as my zeolite box is non-ortho and I have to define 3x3 cell lengths as well
    # box 2 is isotropic liq, so it just needs celllengths (1x3)
    # change your cell vectors and cell sides here
    CellLengths[0][0][0] = lengths1[0] * 10
    CellLengths[0][0][1] = lengths1[1] * 10
    CellLengths[0][0][2] = lengths1[2] * 10
    CellLengths[1][0][0] = lengths2[0] * 10
    CellLengths[1][0][1] = lengths2[1] * 10
    CellLengths[1][0][2] = lengths2[2] * 10
    ##########################################################

    atom_list = []
    molecule_names = []
    nbeads_list = []
    charge_list = []
    for i in range(len(molecules)):
        current_molecule = mb.clone(molecules[i])
        molecule_names.append(current_molecule.name)
        nbeads_list.append(0)
        for particle in current_molecule.particles():
            atom_list.append(particle.name)
            nbeads_list[-1] += 1
            charge_list.append(particle.charge)

    AtomsBox = {}

    AtomsBox[1] = atom_list
    AtomsBox[2] = atom_list

    MoleculesBox = {}
    MoleculesBox[1] = molecule_names
    MoleculesBox[2] = molecule_names
    NBeadsBox = {}
    NBeadsBox[
        1
    ] = nbeads_list  # list that contains number of beads for each molecule type
    NBeadsBox[2] = nbeads_list

    NMoleculesBox = {}
    NMoleculesBox[1] = num_each_moltype_box1  # number of molecules
    NMoleculesBox[2] = num_each_moltype_box2  # number of molecules
    charge_Box = {}
    charge_Box[1] = charge_list
    charge_Box[2] = charge_list

    ####Coordinate file names
    fileBox = {}
    fileBox[1] = xyz_file[0]
    fileBox[2] = xyz_file[1]

    #########config file names
    config_file = {}

    for j in range(NBox):
        if lconfig_file[j]:
            number_of_beads_in_config = sum(
                [a * b for a, b in zip(NBeadsBox[j + 1], NMoleculesBox[j + 1])]
            )
            call(
                "tail -n {} ".format(2 * number_of_beads_in_config)
                + config_file[j + 1]
                + "> temp.xyz",
                shell=True,
            )
            call("""awk '{print>"line-"NR%2}' temp.xyz""", shell=True)
            i = 0
            file_object1 = open("line-1", "r")
            file_object2 = open(fileBox[j + 1], "r")

            next(file_object2)
            next(file_object2)

            a = "\n\n"
            while i < number_of_beads_in_config:
                b = file_object2.readline().split()[0]

                c = file_object1.readline()
                d = b + c

                a = a + d
                i = i + 1
            file_object3 = open("coordBox{}.xyz".format(j + 1), "w")
            file_object3.write(a)
            file_object3.close()
            fileBox[j + 1] = "coordBox{}.xyz".format(j + 1)

    nbeads_diction = {}
    charge_diction = {}

    for i in range(NBox):
        charge_diction[i + 1] = dict(zip(AtomsBox[i + 1], charge_Box[i + 1]))

    for i in range(NBox):
        for j in range(len(MoleculesBox[i + 1])):
            nbeads_diction[MoleculesBox[i + 1][j]] = NBeadsBox[i + 1][j]

    nchain = 0
    for i in range(NBox):
        nchain = nchain + sum(NMoleculesBox[i + 1])

    molecule_list = []
    for i in range(NBox):
        molecule_list.extend(MoleculesBox[i + 1])
    molecule_list = unique(molecule_list)

    nmolty = len(molecule_list)
    nunit = [nbeads_diction[molecule] for molecule in molecule_list]

    liq = {}
    for box in range(1, NBox + 1):
        i = 0  # counter for total atoms
        ni = 0  # counter for atoms in each molecule
        totalbeads = sum(
            [a * b for a, b in zip(NBeadsBox[box], NMoleculesBox[box])]
        )
        if totalbeads == 0:
            continue
        liq[box] = [[99999 for x in range(6)] for y in range(totalbeads)]

        number_type_molecules = NMoleculesBox[box]
        number_beads_each_type = NBeadsBox[box]

        with open(fileBox[box], "r") as f:
            next(f)
            next(f)
            for line in f:
                list1 = line.split()
                element = list1[0]
                x = float(list1[1])
                y = float(list1[2])
                z = float(list1[3])

                liq[box][int(i)][0] = element
                incremental_list = []
                number = 0
                i = i + 0.01
                for j in range(len(number_type_molecules)):
                    number += (
                        number_type_molecules[j] * number_beads_each_type[j]
                    )
                    incremental_list.append(number)
                incremental_list.append(i)
                incremental_list.sort()
                #           print(incremental_list)

                type_element_i = incremental_list.index(i)
                #           print(type_element_i)
                i = int(i)
                liq[box][i][1] = MoleculesBox[box][type_element_i]
                # print(MoleculesBox[box][type_element_i])
                liq[box][i][2] = x
                liq[box][i][3] = y
                liq[box][i][4] = z
                liq[box][i][5] = charge_diction[box][element]
                i = i + 1
        f.close()

    identity_molecule = []
    identity_molecule_box = []
    molecule_number = [i + 1 for i in range(nmolty)]
    molecule_identifier_diction = dict(zip(molecule_list, molecule_number))
    print(molecule_identifier_diction)
    for j in range(NBox):
        if True:
            box = j + 1
            number_type_molecules = NMoleculesBox[box]
            number_beads_each_type = NBeadsBox[box]
            incremental_list = []

            number = 0
            k = 0

            while k < len(number_type_molecules):
                l = 0
                while l < (number_type_molecules[k]):
                    m = 0
                    # while m<((number_beads_each_type[k])):
                    number += number_beads_each_type[k]
                    # m+=1
                    incremental_list.append(int(number))
                    l += 1
                k += 1
            incremental_list = [x - 1 for x in incremental_list]
            for i in incremental_list:
                identity_molecule.append(
                    molecule_identifier_diction[liq[j + 1][i][1]]
                )
                identity_molecule_box.append(j + 1)

    f = open(output_file, "w")
    f.write("           0\n")
    f.write("{0:24.12f}{1:24.12f}{2:24.12f}\n".format(0.01, 0.01, 0.01))
    for i in range(0, NBox):  # number of boxes
        for j in range(nmolty):
            f.write(
                "{0:24.12f}{1:24.12f}{2:24.12f}\n".format(0.3, 0.3, 0.3)
            )  # translation
            f.write(
                "{0:24.12f}{1:24.12f}{2:24.12f}\n".format(0.1, 0.1, 0.1)
            )  # rotation

    a = int(nmolty / 3)
    b = nmolty % 3
    for j in range(NBox):
        for i in range(a):
            f.write("{0:24.12f}{1:24.12f}{2:24.12f}\n".format(5.1, 0.1, 0.1))

        if b == 1:
            f.write(
                "{0:24.12f}\n".format(0.1)
            )  # fluctuating charges for molecule 7
        elif b == 2:
            f.write(
                "{0:24.12f}{1:24.12f}\n".format(0.1, 0.1)
            )  # fluctuating charges for molecule 7

    volume_disp = [1] * NBox
    for element in volume_disp:
        f.write(
            "{0:24.12f}".format(element)
        )  # max volume displacements for boxes 1 - 3
    f.write("\n")

    for j in range(NBox):
        if not lortho[j]:
            a = CellVec[j].tolist()
            for row in a:
                for element in row:
                    f.write("{0:24.12f}".format(element))
                f.write("\n")
            a = CellLengths[j].tolist()

            f.write("{0:24.12f}{1:24.12f}{2:24.12f}\n".format(a[0][0], 0, 0))
            f.write("{0:24.12f}{1:24.12f}{2:24.12f}\n".format(0, a[1][1], 0))
            f.write("{0:24.12f}{1:24.12f}{2:24.12f}\n".format(0, 0, a[2][2]))

        else:
            a = CellLengths[j].tolist()[0]
            f.write("{0:24.12f}{1:24.12f}{2:24.12f}\n".format(a[0], a[1], a[2]))

    f.write("{0:12.0f}\n".format(nchain))  # Total number of molecules
    f.write("{0:12.0f}\n".format(nmolty))  # Total number of types of molecules

    for element in nunit:
        f.write("{0:12.0f}".format(element))
    f.write("\n")

    a = int(nchain / 6)
    b = nchain % 6

    element_index = 0
    for element in identity_molecule:
        f.write("{0:12.0f} ".format(element))
        element_index += 1
        if element_index % 6 == 0 and not (
            element_index == len(identity_molecule)
        ):
            f.write("\n")
    f.write("\n")

    element_index = 0
    for element in identity_molecule_box:
        f.write("{0:12.0f} ".format(element))

        element_index += 1
        if element_index % 6 == 0 and not (
            element_index == len(identity_molecule_box)
        ):
            f.write("\n")
    f.write("\n")

    for j in range(0, NBox):
        for i in range(
            0,
            (
                sum(
                    [
                        a * b
                        for a, b in zip(NBeadsBox[j + 1], NMoleculesBox[j + 1])
                    ]
                )
            ),
        ):
            f.write(
                "{0:24.16f}{1:24.16f}{2:24.16f}\n{3:24.16f}\n".format(
                    liq[j + 1][i][2],
                    liq[j + 1][i][3],
                    liq[j + 1][i][4],
                    liq[j + 1][i][5],
                )
            )

    f.close()

    if os.path.isfile("line-1"):
        os.remove("line-1")

    if os.path.isfile("line-0"):
        os.remove("line-0")

    if os.path.isfile("temp.xyz"):
        os.remove("temp.xyz")

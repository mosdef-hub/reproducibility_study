"""Helper functions for fort77maker_onebox.py and fort77maker_twobox.py."""


def unique(list1):
    """Delete duplicate elements in a list.

    Parameters
    ----------
    list1 : Input list that may contain duplicate elements

    Returns
    -------
    unique_list : Output list with duplicate elements deleted

    """
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
    return unique_list


# This function takes name_list, coordinates, and writes the xyz file in out_filename
def xyzwriter(name_list, coordinates, out_filename):
    """Write an xyz file from element names and coordinates.

    Parameters
    ----------
    name_list : List containing element names
    coordinates : xyz coordinates of elements in AA
    out_filename : Name of the output filename
    """
    if not (len(name_list) == coordinates.shape[0]):
        print(
            "Error in xyzwriter. The atom list size and coordinate size is not same"
        )
        raise IndexError("The atom and coordinate lengths are different.")
    output_string = ""
    output_string += "{}\n".format(len(name_list))
    output_string += "\n"
    for i in range(len(name_list)):
        output_string += "{} {} {} {}\n".format(
            name_list[i],
            coordinates[i, 0],
            coordinates[i, 1],
            coordinates[i, 2],
        )
    with open(out_filename, "w") as text_file:
        text_file.write(output_string)
    print("File written as {}".format(out_filename))

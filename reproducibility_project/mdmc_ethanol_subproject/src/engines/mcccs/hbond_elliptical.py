import os
import freud
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import matplotlib

def getGroups(traj, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail):
    """Get the groups involved in Hydrogen bonding.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    sel_oxygen_head: mdtraj selection string
        String for selecting the acceptor oxygens from the topology file.
    sel_oxygen_tail: mdtraj selection string
        String for selecting the  oxygens bonded to the donor Hydrogen atom from the topology file.
    sel_hydrogen: mdtraj selection string
        String for selecting the donor Hydrogen atoms from the topology file.
    list_names_hydrogen: List of strings
        List of strings containing the names (in the top file) of the Hydrogen atoms making H-bonds
    list_names_oxygen_head: List of strings
        List of strings containing the names (in the top file) of the acceptor Oxygen atoms
    list_names_oxygen_tail: List of strings
        List of strings containing the names (in the top file) of the oxygen atoms bonded to the donor H-atoms

    Returns
    -------

    OxygenHead : List of integres
        List of acceptor oxygen atoms.
    Hydrogen : List of integers
        List of donor Hydrogen atoms.
    OxygenTail : List of integers
        List of oxygen atoms bonded to the donor H atom.
    list_names_hydrogen : List of strings
        List of the names of the donor H atoms.
    list_names_oxygen_head : List of strings
        List of the names of the acceptor O atoms.
    list_names_oxygen_tail : List of strings
        List of the names of the O atoms boned to donor H atoms.
    """



    top = traj.topology
    OxygenHead = top.select(sel_oxygen_head)
    Hydrogen = top.select(sel_hydrogen)
    OxygenTail = top.select(sel_oxygen_tail)

    return OxygenHead, Hydrogen, OxygenTail, list_names_hydrogen,list_names_oxygen_head, list_names_oxygen_tail

def ellipticalFun(dist,degangle):
    """Check if a point falls in the elliptical region.
    Parameters
    ----------
    dist : float
        Distance in nm
    degangle : float
        Angle in degrees

    Returns
    -------
    Boolean
        Returns True if the point lies in the elliptical region. Otherwise False.
    """

    if ((((dist-0.275)/0.05)**2 + (1/50**2)*(degangle-180)**2) < 1): return True

    else: return False

def create_bond_dict(traj, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail, bonded_pdb_provided=False, max_r_bond=0.12):
    """If the bond information is not provided in the traj.top file, then this function helps create bonds between hydrogens and oxygens based on a distance cutoff.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    sel_oxygen_head: mdtraj selection string
        String for selecting the acceptor oxygens from the topology file.
    sel_oxygen_tail: mdtraj selection string
        String for selecting the  oxygens bonded to the donor Hydrogen atom from the topology file.
    sel_hydrogen: mdtraj selection string
        String for selecting the donor Hydrogen atoms from the topology file.
    list_names_hydrogen: List of strings
        List of strings containing the names (in the top file) of the Hydrogen atoms making H-bonds
    list_names_oxygen_head: List of strings
        List of strings containing the names (in the top file) of the acceptor Oxygen atoms
    list_names_oxygen_tail: List of strings
        List of strings containing the names (in the top file) of the oxygen atoms bonded to the donor H-atomsi
    bonded_pdb_provided : Boolean
        Boolean containing the info if a PDB with bond information is provided on not
    max_r_bond : float
        Distance (nm) cutoff for bond between O and H

    Returns
    -------
    OxygenTail_bond_diction : Dictionary
        Dictory containing the information about bond between atoms so that the information is required to be retrieved from the top inforamation in the main loop
    top : mdTraj topology
        Updated topology after bond information updation
    """




    top = traj.top
    atoms = list(top.atoms)
    OxygenHead, Hydrogen, OxygenTail, Hydrogen_atom_names, OxygenHead_names, OxygenTail_names = getGroups(traj,sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
    OxygenTail_bond_diction={}

    if bonded_pdb_provided == False:
        for OxygenTail_index in  OxygenTail:

            OxygenTail_bond_diction[OxygenTail_index]=[]

        frame = 0
        L = traj[frame].unitcell_lengths
        #print(L)
        box = freud.box.Box(Lx=L[0][0],Ly= L[0][1], Lz=L[0][2])
        points_O_tail=traj.xyz[frame][OxygenTail]
        points_H = traj.xyz[frame][Hydrogen]
        aq = freud.locality.AABBQuery(box, points_H)


        query_result = aq.query(points_O_tail, dict(r_max=max_r_bond))
        for bond in query_result:
            if np.isclose(bond[2], 0): #means that the distance is zero
                continue
            point = bond[0]
            neighbor = bond[1]
            dist = bond[2]
            #print(point, neighbor)
            original_Otail_index = OxygenTail[point]
            original_H_index = Hydrogen[neighbor]
            OxygenTail_bond_diction[original_Otail_index].append(original_H_index)
            top.add_bond(atoms[original_Otail_index],atoms[original_H_index ])



        return OxygenTail_bond_diction, top



    for OxygenTail_index in  OxygenTail:

        OxygenTail_bond_diction[OxygenTail_index]=[]
        for bond in top.bonds:

            if (bond.atom1.index == OxygenTail_index) or (bond.atom2.index == OxygenTail_index):

                if bond.atom1.name in Hydrogen_atom_names:
                    bonded_H_index = bond.atom1.index

                elif bond.atom2.name in Hydrogen_atom_names:
                    bonded_H_index = bond.atom2.index

                else:
                    continue
                OxygenTail_bond_diction[OxygenTail_index].append(bonded_H_index)

            else:
                continue
    return OxygenTail_bond_diction, top


def calualateHBMap(traj, r_cutoff, nbins_r, nbins_a, skip_every_x_frames, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail, bonded_pdb_provided = False):
    """Function for calculating the 2D histrogram (cos angle vs distance) for determing H bonds.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    r_cutoff : float
        Distance (nm) up to which pairs are considered.
    nbins_r : int
        Number of bins in the r direction
    nbins_a : int
        Number of bins in the angle direction
    skip_ever_x_frames : int
        Number of frames to be skipped in each loop
    sel_oxygen_head: mdtraj selection string
        String for selecting the acceptor oxygens from the topology file.
    sel_oxygen_tail: mdtraj selection string
        String for selecting the  oxygens bonded to the donor Hydrogen atom from the topology file.
    sel_hydrogen: mdtraj selection string
        String for selecting the donor Hydrogen atoms from the topology file.
    list_names_hydrogen: List of strings
        List of strings containing the names (in the top file) of the Hydrogen atoms making H-bonds
    list_names_oxygen_head: List of strings
        List of strings containing the names (in the top file) of the acceptor Oxygen atoms
    list_names_oxygen_tail: List of strings
        List of strings containing the names (in the top file) of the oxygen atoms bonded to the donor H-atoms
    bonded_pdb_provided : Boolean
        Boolean indicating if the topology file contains the bond information

    Returns
    -------
    rdf_output : numpy ndarray
        bins for r (nm), g(r)
    inter_output : numpy ndarray
        bins for angle, adf
    map_output : numpy 2D array
        2D array for the frequency on the 2D grid
    ans_overal : list of nhbond_i x 3 ndarray with a length equal to numframes 
        nhbonds_i x 3 array that contains the pair of indices of atoms that make a H bond (O_Tail, H, O_Head)
    hbond_time : List of length equal to numframes
        number of Hbonds in each frame


    Examples
    --------
    >>> nbins_r = 200 ; nbins_a = 200 ; r_cutoff = 0.75 ; skip_every_x_frames = 1
    >>> sel_oxygen_head = 'name O5' ; sel_hydrogen = 'name H1 or name H2' ; sel_oxygen_tail = 'name O5'; list_names_hydrogen = ["H1", "H2"] ; list_names_oxygen_head = ["O5"] ; list_names_oxygen_tail = ["O5"]
    >>> rdf_output, inter_output, map_output,hbond,hbond_time = calualateHBMap(traj, r_cutoff, nbins_r, nbins_a, skip_every_x_frames, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
    >>> #Plotting
    >>> plt.figure() ; cmap = plt.get_cmap('jet') ; plt.figure(figsize=(5, 3)) ; plt.style.use('default'); levels = np.linspace(0,10,11) ; cs = plt.contourf(rdf_output[0], inter_output[0], map_output,levels = levels, cmap=cmap) ; plt.xlabel('r (nm)') ; plt.ylabel('\u03B8 (degrees)') ; plt.xlim([0.2, 0.4]) ; plt.ylim([140, 180]) ; plt.colorbar()

    """

    bond_diction, top = create_bond_dict(traj, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail, bonded_pdb_provided = bonded_pdb_provided)
    traj.top = top
    top = traj.top
    print("traj as {} bonds".format(top.n_bonds))
    radii, dr = np.linspace(0, r_cutoff, nbins_r+1, retstep=True)

    angles, da = np.linspace(0, np.pi, nbins_a+1, retstep=True)

    print("da", da)
    print("dr", dr)
    volume_shells = np.zeros(nbins_r);
    for i in range(0, nbins_r):
        volume_shells[i] = (4/3)*np.pi*(radii[i+1]**3-radii[i]**3)

    r_centers = 0.5 * (radii[1:] + radii[:-1])
    a_centers = 0.5 * (angles[1:] + angles[:-1])

    """ Initialize output array """
    gr_final = np.zeros(nbins_r)
    intangle_final = np.zeros(nbins_a)
    maps_final = np.zeros((nbins_a,nbins_r))
    hbond_time = []
    ans_overall = []
    prev = 0
    frames_calculated = 0
    """ Iterate through each frame """
    for frame in range(0, traj.n_frames, skip_every_x_frames):
        local_frame = traj[frame]
        ans = np.array([]).reshape(0,3)
        frames_calculated +=1
        print("working on {}".format(frame))

        L = traj[frame].unitcell_lengths
        #print(L)
        box = freud.box.Box(Lx=L[0][0],Ly= L[0][1], Lz=L[0][2])
        intangle_prob = np.zeros(nbins_a)
        g_of_r = np.zeros(nbins_r)
        maps = np.zeros((nbins_a, nbins_r))

        """ Get Neighbors """
        OxygenHead, Hydrogen, OxygenTail, Hydrogen_atom_names, OxygenHead_names, OxygenTail_names = getGroups(traj,  sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
        numHead = len(OxygenHead)
        numTail = len(OxygenTail)
        #print(numHead, numTail)
        points_O_tail=traj.xyz[frame][OxygenTail]
        points_O_head=traj.xyz[frame][OxygenHead]
        aq = freud.locality.AABBQuery(box, points_O_head)



        """ Iterate through Neighbors """

        query_result = aq.query(points_O_tail, dict(r_max=r_cutoff))
        triplets_analyzed = []
        #n_pairs = 0
        for shbond in query_result:
            if np.isclose(shbond[2], 0): #means that the distance is zero
                continue
            point = shbond[0]
            neighbor = shbond[1]
            dist = shbond[2]
            #print(point, neighbor)
            original_Otail_index = OxygenTail[point]
            original_Ohead_index = OxygenHead[neighbor]
            potential_lists = []
            for bonded_H in bond_diction[original_Otail_index]:
                potential_lists.append([original_Otail_index, bonded_H, original_Ohead_index])
            triplet_already_analyzed = False
            for triplet in potential_lists:
                if triplet in triplets_analyzed:
                    triplet_already_analyzed = True

            if triplet_already_analyzed:
                continue
            else:
                for bonded_H in bond_diction[original_Otail_index]:

                    triplets_analyzed.append([original_Otail_index, bonded_H, original_Ohead_index])
                    triplets_analyzed.append([original_Ohead_index, bonded_H, original_Otail_index])
            #now check if tail is connected to a H
            #print(bond_diction)
            for bonded_H in bond_diction[original_Otail_index]:

                original_H_index = bonded_H


               # print("{}-{} is bonded to {}-{}".format(original_Otail_index,bond.atom1.name, original_H_index,bond.atom2.name))
                index_d = np.digitize(dist, radii)-1
                g_of_r[index_d] += 1.00
                """ calculate angle """
                normal_angle = md.compute_angles(local_frame, np.array([[original_Otail_index, original_H_index, original_Ohead_index]]))[0][0]
                if normal_angle > np.pi:

                    print("Angle is {} radians greater than {}", normal_angle-np.pi, np.pi)
                    normal_angle = normal_angle -1e-3

                index_i = np.digitize(normal_angle,angles)-1
                #print("index_i is",index_i)
                intangle_prob[index_i] += 1.00
                """ put into map """
                maps[index_i,index_d] += 1.00
                deg_angle = 180*normal_angle/np.pi
                #print("{}-{}-{} has a distance of {} and an angle of {}".format(original_Otail_index, original_H_index, original_Ohead_index, dist, deg_angle))
                if ellipticalFun(dist, deg_angle) :
                    #print(np.array([original_Otail_index, original_H_index, original_Ohead_index]), normal_angle)
                    ans = np.vstack([ans, np.array([original_Otail_index, original_H_index, original_Ohead_index])])
                    


        hbond_time.append(ans.shape[0])
        ans_overall.append(ans)
        """ normalization """
        #print("The num pairs found in this frame is {}".format(n_pairs))
        NbyV = numHead * (numTail-1) / (L[0][0]*L[0][1]*L[0][2])
        #print("NbyV is ", NbyV)
        #total_g = numHead * (numTail) / (L[0][0]*L[0][1]*L[0][2])
        total_m = np.zeros(nbins_a)

        """ map normalization """
        for i, value in enumerate(g_of_r):
            g_of_r[i] = value /( 4*np.pi*(r_centers[i]**2)*dr *NbyV)

        for i in range(nbins_a):
            intangle_prob[i] = intangle_prob[i]/np.sin(a_centers[i])

        for i in range(0, nbins_a):
            total_m[i] = sum(maps[i,:])
            if total_m[i] == 0:
                maps[i,:] = 0
                continue

            for j in range(0, nbins_r):
                rho_r_theta = 4*np.pi*(r_centers[j]**2)*abs(np.sin(a_centers[i]))*dr*da

#                 if maps[i,j] > 0:
#                     print("r is ",r_centers[j] )
#                     print("rho_r_theta is ", rho_r_theta)
#                     print("element volume is ", dr*da)
#                     print("freq of this square is", maps[i,j])
#                     print("sin theta is ",np.sin(a_centers[i]) )
#                     print("updated maps[i,j]", maps[i,j] / (rho_r_theta*NbyV))
                maps[i,j] = maps[i,j] / (rho_r_theta*NbyV)



        #g_of_r /= sum(g_of_r)*dr
        intangle_prob = intangle_prob/(da*sum(intangle_prob))

        intangle_final += intangle_prob
        maps_final += maps
        gr_final += g_of_r

    intangle_final /= frames_calculated
    maps_final /= frames_calculated
    gr_final /= frames_calculated

    inter_output = np.array([a_centers, intangle_final])
    rdf_output = np.array([r_centers, gr_final])
    map_output = maps_final
    return rdf_output, inter_output, map_output,ans_overall,hbond_time



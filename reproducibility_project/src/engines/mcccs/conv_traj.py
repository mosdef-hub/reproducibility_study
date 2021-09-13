import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
import shutil
from scipy import stats
import freud
from glob import glob

def main():
    data_path = "analysis_data"
    #if os.path.exists(data_path):
    #    shutil.rmtree(data_path)
    #os.makedirs(data_path)
    #delete the folder manually
    if not os.path.isdir(data_path):

        os.makedirs(data_path)

    os.chdir(data_path)
    

    project = signac.get_project()

    for (molecule, ensemble, temperature, pressure), group in project.groupby(("molecule","ensemble", "temperature","pressure", )):
        print(molecule, ensemble, temperature, pressure)
        if not os.path.isdir("{}_{}_{}K_{}kPa".format(molecule, ensemble, temperature, pressure)):

            os.makedirs("{}_{}_{}K_{}kPa".format(molecule, ensemble, temperature, pressure))
        os.chdir("{}_{}_{}K_{}kPa".format(molecule, ensemble, temperature, pressure))
        base_dir=os.getcwd() 
        if ensemble == "NPT" and molecule == "methaneUA":
            if os.path.exists(os.path.join(os.getcwd(), 'comb_traj.h5')): 
                continue
            traj_list=[]
            for job in group:
                 #print(job)
                 os.chdir(job.ws)
                 traj_files = sorted(glob("box1movie1a*prod*pdb*"))
                 if len(traj_files)<4:
                     #print("warning, only {}  prod cycles complete for {} {} {} {} {}".format(len(traj_files),job, molecule, ensemble, temperature, pressure))
                     print("warning, only {}  prod cycles complete for {} {} {} {} {}. So not inlcuding its xyz files.".format(len(traj_files),job, molecule, ensemble, temperature, pressure))
                     continue
                 prod_number=0
                 for filename in traj_files:
                     print("The filename is {} and the prod number is {}. These two should match.".format(filename,prod_number))
                     traj = md.load(filename, top='init1.mol2')
                     fort12_filename="fort.12.prod{}".format(prod_number)
                     fort12=np.genfromtxt(fort12_filename, skip_header=1)
                     traj= md.Trajectory(traj.xyz, traj.top,unitcell_lengths = fort12[9::10,0:3], unitcell_angles = np.tile([90.,90.,90.], (traj.n_frames,1)),)
                     traj_list.append(traj)
                     print("one traj loaded")
                     prod_number+=1
            
            if len(traj_list)!=64:
                continue
            comb_traj=md.join(traj_list)
            print("A total of {} trajs. Ideal number is be 64.".format(len(traj_list)))
            os.chdir(base_dir)
            comb_traj.save_hdf5('comb_traj.h5')
            
        elif ensemble == "GEMC-NVT" and molecule == "methaneUA":
            for job in group:
                print(job)
                
        os.chdir("..")

if __name__ == "__main__":
    main()

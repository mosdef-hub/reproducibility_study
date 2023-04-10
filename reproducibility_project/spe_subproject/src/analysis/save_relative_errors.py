import signac
import numpy as np
import pandas as pd
p = signac.get_project("./")
df = pd.DataFrame(columns=["jobid", "molecule", "engine", "temperature", "energy"])
for job in p:
    if job.sp.engine == "lammps-UD":
        continue
    if job.sp.engine == "hoomd":
        jobdata = np.genfromtxt(job.ws+"/log-spe-raw.txt", delimiter=" ", dtype=float, skip_header=1)
        pe = jobdata[-1]
    elif job.sp.engine == "gomc":
        jobdata = np.genfromtxt(job.ws+"/log-spe.txt", delimiter=",", dtype=float, skip_header=1)
        pe = jobdata[1]
    else:
        jobdata = np.genfromtxt(job.ws+"/log-spe.txt", delimiter=",", dtype=float, skip_header=1)
        pe = jobdata[0]

    df.loc[len(df.index)] = [job.id, job.sp.molecule, job.sp.engine, job.sp.temperature, pe]

re_func = lambda x: (x-x.mean())/x.mean()*100000
df["relativeError"] = df.groupby(["molecule", "temperature"])["energy"].transform(func=re_func)

df.to_csv("summary_single_point_energies.csv")

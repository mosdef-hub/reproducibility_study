"""Generate CSV from job docs for job densities."""

import os

import camelot
import numpy as np
import pandas as pd
import signac

updated_masses = {
    "methaneUA": 16.043,
    "pentaneUA-flexible_bonds": 72.151,
    "pentaneUA-constrain_bonds": 72.151,
    "pentaneUA": 72.151,
    "benzeneUA": 78.114,
    "waterSPCE": 18.015,  # 18.015324,
    "ethanolAA": 46.069,  # 46.068672,
}


def load_mosdef_jobs():
    """Load MoSDeF densities from job docs."""
    # Load by molecule, temperature, forcefield
    densitiesDict = {}
    p = signac.get_project("./")
    sampled_jobsDict = {}
    for i, job in enumerate(p):
        if not job.isfile("signac_statepoint.json"):
            continue
        filter_list = [
            (
                job.sp.engine == "mcccs"
                and job.sp.mass != updated_masses[job.sp.molecule]
            ),
            job.sp.engine == "lammps-UD",
            job.sp.ensemble != "NPT",
            job.sp.molecule == "pentaneUA",
        ]
        if np.any(filter_list):
            continue

        print(f"{i/len(p)*100:.1f}% done")
        try:
            sampled_jobsDict[job] = (
                job.doc["npt_density_avg"],
                job.doc["npt_density_std"],
            )
            # job.doc['npt/sampling_results']["density"]
        except:
            return job.id, job.doc, job.sp
    return sampled_jobsDict


headers = [
    "associated_work",
    "MCorMD",
    "engine",
    "molecule",
    "replicate",
    "temperature",
    "density",
    "density-std",
    "forcefield",
]
MCMDDict = dict()
for engines, enginetype in list(
    zip(
        [("cassandra", "gomc", "mcccs"), ("lammps-VU", "gromacs", "hoomd")],
        ("MC", "MD"),
    )
):
    for engine in engines:
        MCMDDict.update({engine: enginetype})


def generate_mosdef_job_csv():
    """Create a csv of all data for mosdef jobs based on statepoint to read into dataframes."""
    if os.path.exists("csvs/job_density_data.csv"):
        print("Density data is already generated at csvs/job_density_data.csv.")
        return
    job_densityDict = load_mosdef_jobs()
    row_list = []
    for job in job_densityDict:
        inputs = [
            "MoSDeF",
            MCMDDict[job.sp.engine],
            job.sp.engine,
            job.sp.molecule,
            job.sp.replica,
            job.sp.temperature,
            job_densityDict[job][0],
            job_densityDict[job][1],
            job.sp.forcefield_name,
        ]
        row_list.append(dict((key, val) for key, val in zip(headers, inputs)))
    df = pd.DataFrame(row_list)
    df.to_csv("csvs/job_density_data.csv")


def _generate_csv_from_PDF(molecule):
    """Methods to generate dataframe from PDF supplied for paper SI."""
    page = molecule_pageDictionary[molecule]

    # Read remote pdf into a list of DataFrame
    link = "Hasse_SI.pdf"  # must have access locally for this python file to load it.
    tables = camelot.read_pdf(link, flavor="stream", pages=page, edge_tol=500)
    if len(tables) == 1:
        table = tables[0]
        molecule_nameStr = table.df[1][0]
        if molecule in [
            "OPLS Ethane at 5 MPa",
            "OPLS Propane at 5 MPa",
            "OPLS n-Butane at 5 MPa",
            "OPLS iso-Butane at 5 MPa",
        ]:
            df = table.df[table.df.columns[1:]][2:]
            df.columns = np.arange(len(df.columns))
        elif "AMBER" in molecule:
            df = table.df[table.df.columns[:]][2:]
        else:
            df = table.df[table.df.columns[:]][1:]
        row_start = 1 + np.where(df[0] == "Program(Group)")[0][0]
        headers = df.iloc[row_start - 1]
        newDF = pd.DataFrame(df.values[row_start:], columns=headers)
    else:
        # need to combine table0 and table1 parts
        df = pd.concat((tables[0].df, tables[1].df), axis=0)
        header_ids = np.where(df[0] == "Program(Group)")[0]
        if len(header_ids) == 2:
            row_start = 1 + header_ids[0]
            row_end = header_ids[1] - 2
            headers = df.iloc[row_start - 1]
            newDF = pd.DataFrame(df.values[row_start:row_end], columns=headers)
        elif len(header_ids) == 1:
            row_start = 1 + header_ids[0]
            headers = df.iloc[row_start - 1]
            newDF = pd.DataFrame(df.values[row_start:], columns=headers)

    for i in range(len(newDF.index)):
        headerInt = int(4 * np.floor(i / 4))
        header = newDF.iloc[headerInt]["Program(Group)"]
        newDF.iloc[i]["Program(Group)"] = header

    # we will change the data type
    # of id column to str by giving
    # the dict to the astype method
    columns = list(newDF.columns[1:])

    for c in columns:
        newDF[c] = pd.to_numeric(newDF[c], errors="coerce")
    newDF["statepoint"] = "temp: " + newDF["T / K"].astype(str)
    newDF["Relative_Error"] = newDF.groupby("statepoint")[
        "ρ / kg m−3"
    ].transform(_calculate_relative_error)
    newDF.to_csv(f"csvs/{molecule}_data.csv")


def _calculate_relative_error(ser):
    """Calculate relative density error for all engines at that statepoint."""
    dataArray = ser.to_numpy()
    mean = np.nanmean(dataArray)
    return (dataArray - mean) / mean * 100


if __name__ == "__main__":
    generate_mosdef_job_csv()

    molecule_pageDictionary = {
        "OPLS Ethane at 5 MPa": "8",
        "OPLS Ethane at 41 MPa": "9",
        "OPLS Ethane at 70 MPa": "10",
        "TraPPE Ethane at 5 MPa": "11",
        "TraPPE Ethane at 41 MPa": "12",
        "TraPPE Ethane at 70 MPa": "13",
        "OPLSAMBER Ethane at 5 MPa": "14,15",
        "OPLSAMBER Ethane at 41 MPa": "15,16",
        "OPLSAMBER Ethane at 70 MPa": "16,17",
        "OPLS Propane at 5 MPa": "18",
        "OPLS Propane at 41 MPa": "19",
        "OPLS Propane at 70 MPa": "20",
        "TraPPE Propane at 5 MPa": "21",
        "TraPPE Propane at 41 MPa": "22",
        "TraPPE Propane at 70 MPa": "23",
        "OPLSAMBER Propane at 5 MPa": "24-25",
        "OPLSAMBER Propane at 41 MPa": "25-26",
        "OPLSAMBER Propane at 70 MPa": "26-27",
        "OPLS n-Butane at 5 MPa": "28",
        "OPLS n-Butane at 41 MPa": "29",
        "OPLS n-Butane at 70 MPa": "30",
        "TraPPE n-Butane at 5 MPa": "31",
        "TraPPE n-Butane at 41 MPa": "32",
        "TraPPE n-Butane at 70 MPa": "33",
        "OPLSAMBER n-Butane at 5 MPa": "34,35",
        "OPLSAMBER n-Butane at 41 MPa": "35,36",
        "OPLSAMBER n-Butane at 70 MPa": "36,37",
        "OPLS iso-Butane at 5 MPa": "38",
        "OPLS iso-Butane at 41 MPa": "39",
        "OPLS iso-Butane at 70 MPa": "40",
        "TraPPE iso-Butane at 5 MPa": "41",
        "TraPPE iso-Butane at 41 MPa": "42",
        "TraPPE iso-Butane at 70 MPa": "43",
        "OPLSAMBER iso-Butane at 5 MPa": "44,45",
        "OPLSAMBER iso-Butane at 41 MPa": "45,46",
        "OPLSAMBER iso-Butane at 70 MPa": "46,47",
    }
    for molecule in molecule_pageDictionary:
        print(f"Generating CSV for: {molecule}")
        print("################################")
        print("################################")
        _generate_csv_from_PDF(molecule)
    print("\n\nCompleted Generating all CSVs in the directory csvs/\n\n")

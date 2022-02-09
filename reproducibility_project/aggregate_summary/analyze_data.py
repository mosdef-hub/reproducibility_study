"""Create summary data based on the aggregate values."""
import pandas as pd
import signac


def prepare_data(p: signac.Project) -> pd.DataFrame:
    """Prepare the data in the aggregate summary, return a pandas dataframe."""
    summary = []
    header = [
        "molecule",
        "engine",
        "temperature",
        "pressure",
        "ensemble",
        "forcefield_name",
        "cutoff_style",
        "long_range_correction",
        "r_cut",
        "pressure-avg",
        "pressure-std",
        "pressure-sem",
        "temperature-avg",
        "temperature-std",
        "temperature-sem",
        "potential_energy-avg",
        "potential_energy-std",
        "potential_energy-sem",
        "density-avg",
        "density-std",
        "density-sem",
    ]

    j_list = []
    for j in p:
        if j.doc:
            j_list.append(j)

    for j in j_list:
        tmp = []
        for col in header:
            if col in j.sp:
                tmp.append(j.sp[col])
            elif col in j.doc:
                tmp.append(j.doc[col])
            else:
                tmp.append(None)
        summary.append(tmp)
        del tmp

    df = pd.DataFrame(summary, columns=header)
    return df


def main():
    """Generate summary files for the aggregate project."""
    p = signac.get_project()

    molecule_set = set()
    for job in p:
        molecule_set.add(job.sp["molecule"])

    summary_df = prepare_data(p)
    summary_df.to_csv("aggregate_summary_all.csv")

    group_key = "molecule"
    for molecule in molecule_set:
        fname = f"aggregate_summary_{molecule}.csv"
        try:
            mol_group = summary_df.groupby(group_key)
            mol_df = mol_group.get_group(molecule)
        except KeyError:
            print(f"skipping: {molecule}, no data available.")
            continue
        mol_df.to_csv(fname)


if __name__ == "__main__":
    main()

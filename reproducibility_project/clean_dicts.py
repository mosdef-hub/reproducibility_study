"""Remove job doc entries related to the analysis of the simulation data."""

import signac


def main():
    """Remove entries from the job docs based on specific keys."""
    proj = signac.get_project()
    print(proj.id)

    keys_to_pop = [
        "npt/max_t0",
        "nvt/max_t0",
        "nvt/sampling_results",
        "npt/sampling_results",
    ]
    prefixes = ["npt", "nvt"]
    suffixes = ["avg", "std"]
    props = [
        "temperature",
        "potential_energy",
        "kinetic_energy",
        "pressure",
        "density",
        "volume",
    ]

    for prefix in prefixes:
        for suffix in suffixes:
            for prop in props:
                keys_to_pop.append(f"{prefix}_{prop}_{suffix}")

    print(keys_to_pop)

    for j in proj:
        for key in keys_to_pop:
            if j.doc.get(key) is not None:
                _ = j.doc.pop(key)


if __name__ == "__main__":
    main()

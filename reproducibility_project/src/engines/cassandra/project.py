"""Setup for signac, signac-flow, signac-dashboard for this study."""
# import foyer
import os
from pathlib import Path

import flow
from flow import directives

import reproducibility_project.templates.ndcrc


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


@Project.label
def is_cassandra(job):
    """Verify that the job is for Cassandra."""
    return job.sp.engine == "cassandra"


@Project.label
def cassandra_complete(job):
    """Check whether job is complete."""
    return check_complete(job.fn("prod"))


def check_complete(run_name):
    """Check whether MoSDeF Cassandra simulation with run_name has completed."""
    complete = False
    fname = run_name + ".out.log"
    if not os.path.exists(fname):
        return complete
    with open(fname) as f:
        for line in f:
            if "Cassandra simulation complete" in line:
                complete = True
                break
    return complete


@Project.operation
@Project.pre(is_cassandra)
@Project.post(cassandra_complete)
@directives(omp_num_threads=4)
def run_cassandra(job):
    """Run all simulation stages for given statepoint."""
    import foyer
    import mbuild as mb
    import mosdef_cassandra as mc
    import numpy as np
    import unyt as u
    from mbuild.formats.xyz import read_xyz
    from pymbar.timeseries import detectEquilibration

    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
        get_molecule,
    )
    from reproducibility_project.src.utils.forcefields import load_ff

    molecule = job.sp.molecule

    compound = get_molecule(job.sp)

    ensemble = job.sp.ensemble

    cass_ensembles = {
        "NPT": "npt",
        "GEMC-NVT": "gemc",
    }

    Nliq = job.sp.N_liquid
    if ensemble == "GEMC-NVT":
        Nvap = job.sp.N_vap
    else:
        Nvap = 0

    N = Nliq + Nvap

    Nlist = [[Nliq]]

    equil_length = 40000
    prod_length = 120000

    filled_boxes = construct_system(job.sp)

    liqbox_filled = filled_boxes[0]
    if ensemble == "GEMC-NVT":
        vapbox_filled = filled_boxes[1]
        Nlist.append([Nvap])

    ffname = job.sp.forcefield_name
    ff = load_ff(ffname)
    structure = ff.apply(compound)

    Tmelt = 1000.0 * u.K

    T = job.sp.temperature * u.K
    P = job.sp.pressure * u.kPa

    species_list = [structure]
    cutoff = job.sp.r_cut * u.nm

    seedslist = [
        [7860904, 98601355],
        [955793508, 417823039],
        [254420642, 2130720464],
        [58120272, 1465850411],
        [757616664, 1940492980],
        [966844679, 1326087693],
        [1992175335, 1840317929],
        [650725500, 646331893],
        [1204247127, 1521385831],
        [2000946981, 1969870819],
        [1488434295, 648017520],
        [1128424221, 1140005446],
        [58870203, 2133009902],
        [2024564019, 2014788324],
        [133927152, 2052536489],
        [23375750, 1951798462],
    ]

    seeds = seedslist[job.sp.replica]

    cbmc_n_ins = 12
    cbmc_n_dihed = 50

    prop_freq = 10
    coord_freq = 10
    if ff.combining_rule == "geometric":
        comb_rule = "geometric"
    else:
        comb_rule = "lb"

    proplist = [
        "energy_total",
        "volume",
        "nmols",
        "pressure",
        "density",
        "mass_density",
        "energy_angle",
        "energy_dihedral",
        "energy_intravdw",
        "energy_intraq",
        "energy_inter",
        "energy_intervdw",
        "energy_lrc",
        "energy_interq",
        "energy_recip",
        "energy_self",
        "enthalpy",
    ]

    meltsystem_liq = mc.System([liqbox_filled], species_list, [[Nliq]])

    p_volume = 0.01
    p_translate = 0.0
    p_rotate = 0.0
    p_regrow = 0.0
    p_swap = 0.0

    if molecule == "methaneUA":
        if ensemble == "NPT":
            p_translate = 0.99
        else:
            p_translate = 0.7
            p_swap = 0.29
    elif molecule == "pentaneUA" and ensemble == "GEMC-NVT":
        p_swap = 0.2
        p_translate = 0.27
        p_rotate = 0.26
        p_regrow = 0.26
    elif molecule == "waterSPCE":
        p_translate = 0.49
        p_rotate = 0.5
    else:
        p_translate = 0.33
        p_rotate = 0.33
        p_regrow = 0.33

    nvtmoves = mc.MoveSet("nvt", species_list)
    nvtmoves.prob_rotate = p_rotate
    nvtmoves.prob_translate = p_translate
    nvtmoves.prob_regrow = p_regrow

    moveset = mc.MoveSet(cass_ensembles[ensemble], species_list)
    moveset.prob_volume = p_volume
    moveset.prob_translate = p_translate
    moveset.prob_rotate = p_rotate
    moveset.prob_regrow = p_regrow
    moveset.cbmc_n_ins = cbmc_n_ins
    moveset.cbmc_n_dihed = cbmc_n_dihed

    with job:
        if not check_complete("nvt_melt"):
            mc.run(
                system=meltsystem_liq,
                moveset=nvtmoves,
                run_type="equilibration",
                run_length=5000,
                temperature=Tmelt,
                properties=proplist,
                cutoff_style="cut",
                vdw_cutoff=cutoff,
                charge_cutoff=cutoff,
                run_name="nvt_melt",
                prop_freq=prop_freq,
                coord_freq=5000,
                units="sweeps",
                steps_per_sweep=Nliq,
                seeds=seeds,
            )

        nvtendbox_liq = read_xyz("nvt_melt.out.xyz")

        nvtendbox_liq.box = liqbox_filled.box

        nvtsystem_liq = mc.System([nvtendbox_liq], species_list, [[Nliq]])

        if not check_complete("nvt_equil"):
            mc.run(
                system=nvtsystem_liq,
                moveset=nvtmoves,
                run_type="equilibration",
                run_length=5000,
                temperature=T,
                properties=proplist,
                cutoff_style="cut",
                vdw_cutoff=cutoff,
                charge_cutoff=cutoff,
                run_name="nvt_equil",
                prop_freq=prop_freq,
                coord_freq=5000,
                units="sweeps",
                steps_per_sweep=Nliq,
                seeds=seeds,
            )

        nvtendbox_liq = read_xyz("nvt_equil.out.xyz")

        nvtendbox_liq.box = liqbox_filled.box
        boxlist = [nvtendbox_liq]

        if ensemble == "GEMC-NVT":
            meltsystem_vap = mc.System([vapbox_filled], species_list, [[Nvap]])
            if not check_complete("nvt_melt_vap"):
                mc.run(
                    system=meltsystem_vap,
                    moveset=nvtmoves,
                    run_type="equilibration",
                    run_length=5000,
                    temperature=Tmelt,
                    properties=proplist,
                    cutoff_style="cut",
                    vdw_cutoff=cutoff,
                    charge_cutoff=cutoff,
                    run_name="nvt_melt_vap",
                    prop_freq=prop_freq,
                    coord_freq=5000,
                    units="sweeps",
                    steps_per_sweep=Nvap,
                    seeds=seeds,
                )
            meltendbox_vap = read_xyz("nvt_melt_vap.out.xyz")

            meltendbox_vap.box = vapbox_filled.box
            nvtsystem_vap = mc.System([meltendbox_vap], species_list, [[Nvap]])
            if not check_complete("nvt_equil_vap"):
                mc.run(
                    system=nvtsystem_vap,
                    moveset=nvtmoves,
                    run_type="equilibration",
                    run_length=5000,
                    temperature=T,
                    properties=proplist,
                    cutoff_style="cut",
                    vdw_cutoff=cutoff,
                    charge_cutoff=cutoff,
                    run_name="nvt_equil_vap",
                    prop_freq=prop_freq,
                    coord_freq=5000,
                    units="sweeps",
                    steps_per_sweep=Nvap,
                    seeds=seeds,
                )
            nvtendbox_vap = read_xyz("nvt_equil_vap.out.xyz")

            nvtendbox_vap.box = vapbox_filled.box
            boxlist.append(nvtendbox_vap)
            moveset.prob_swap = p_swap

        system = mc.System(boxlist, species_list, Nlist)

        l_equil = False
        cycles_done = 0

        if not check_complete("equil_0"):
            mc.run(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=equil_length,
                temperature=T,
                pressure=P,  # this line is ignored if ensemble isn't NPT
                properties=proplist,
                cutoff_style="cut",
                vdw_cutoff=cutoff,
                charge_cutoff=cutoff,
                run_name="equil_0",
                prop_freq=prop_freq,
                coord_freq=coord_freq,
                units="sweeps",
                steps_per_sweep=N,
                seeds=seeds,
            )

        prior_run = "equil_" + str(cycles_done)
        cycles_done += equil_length

        if ensemble == "GEMC-NVT":
            prpsuffix_liq = ".out.box1.prp"
            prpsuffix_vap = ".out.box2.prp"
        else:
            prpsuffix_liq = ".out.prp"

        t, g, Neff_max = detectEquilibration(
            np.loadtxt(prior_run + prpsuffix_liq, usecols=5)
        )

        if ensemble == "GEMC-NVT":
            tvap, gvap, Neff_max_vap = detectEquilibration(
                np.loadtxt(prior_run + prpsuffix_vap, usecols=5)
            )
            t = max(t, tvap)

        for i in range(2):
            if t >= equil_length * 3 / (prop_freq * 4):
                if not check_complete("equil_" + str(cycles_done)):
                    mc.restart(
                        restart_from=prior_run,
                        run_name="equil_" + str(cycles_done),
                        run_type="equilibration",
                        total_run_length=cycles_done + equil_length,
                    )
                prior_run = "equil_" + str(cycles_done)
                cycles_done += equil_length
                t, g, Neff_max = detectEquilibration(
                    np.loadtxt(prior_run + prpsuffix_liq, usecols=5)
                )

                if ensemble == "GEMC-NVT":
                    tvap, gvap, Neff_max_vap = detectEquilibration(
                        np.loadtxt(prior_run + prpsuffix_vap, usecols=5)
                    )
                    t = max(t, tvap)

                for rmtgt in list(Path(".").glob("equil_*.[xH]*")):
                    os.remove(rmtgt)

            else:
                break

        mc.restart(
            restart_from=prior_run,
            run_name="prod",
            run_type="production",
            total_run_length=cycles_done + prod_length,
        )


if __name__ == "__main__":
    pr = Project()
    pr.main()

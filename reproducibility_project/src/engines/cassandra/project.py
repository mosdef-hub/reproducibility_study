"""Setup for signac, signac-flow, signac-dashboard for this study."""
# import foyer
import os
from pathlib import Path

import flow
import numpy as np
import pandas as pd
from flow import directives

import reproducibility_project.templates.ndcrc


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = Path(os.getcwd()).absolute()
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
    """Check whether MoSDeF Cassandra simulation with run_name or its last restart has completed."""
    complete = False
    fname = run_name + ".out.log"
    loglist = list_with_restarts(fname)
    if not loglist:
        return complete
    with loglist[-1].open() as f:
        for line in f:
            if "Cassandra simulation complete" in line:
                complete = True
                break
    return complete


def list_with_restarts(fpath):
    """List fpath and its restart versions in order as pathlib Path objects."""
    fpath = Path(fpath)
    if not fpath.exists():
        return []
    parent = fpath.parent
    fname = fpath.name
    fnamesplit = fname.split(".out.")
    run_name = fnamesplit[0]
    suffix = fnamesplit[1]
    restarts = [
        Path(parent, f)
        for f in sorted(list(parent.glob(run_name + ".rst.*.out." + suffix)))
    ]
    restarts.insert(0, fpath)  # prepend fpath to list of restarts
    return restarts


def get_last_checkpoint(run_name):
    """Get name of last restart based on run_name."""
    fname = run_name + ".out.chk"
    return list_with_restarts(fname)[-1].name.split(".out.")[0]


def has_checkpoint(run_name):
    """Check whether there is a checkpoint for run_name."""
    fname = run_name + ".out.chk"
    return os.path.exists(fname)


def merge_restart_prp(fname):
    """Merge restart prp files."""
    pathlist = list_with_restarts(fname)
    if len(pathlist) < 2:
        return
    lastcycles = [np.loadtxt(f, usecols=0, dtype=int)[-1] for f in pathlist]
    with pathlist[0].open("a") as mainfile:
        for i in range(1, len(pathlist)):
            with pathlist[i].open() as rfile:
                for j in range(3):
                    next(rfile)
                line = rfile.readline()
                n = int(line.split()[0])
                while n <= lastcycles[i - 1]:
                    line = rfile.readline()
                    n = int(line.split()[0])
                mainfile.write(line)
                for line in rfile:
                    mainfile.write(line)
            os.remove(pathlist[i])


def merge_restart_traj(fname):
    """Merge restart .H or .xyz files."""
    pathlist = list_with_restarts(fname)
    if len(pathlist) < 2:
        return
    with pathlist[0].open("a") as mainfile:
        for i in range(1, len(pathlist)):
            with pathlist[i].open() as rfile:
                for line in rfile:
                    mainfile.write(line)
            os.remove(pathlist[i])


@Project.operation
@Project.pre(is_cassandra)
@Project.post(cassandra_complete)
@directives(omp_num_threads=4)
def run_cassandra(job):
    """Run all simulation stages for given statepoint."""
    import foyer
    import mbuild as mb
    import mosdef_cassandra as mc
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

    cass_cutoffs = {
        ("hard", "None"): "cut",
        ("hard", "energy_pressure"): "cut_tail",
        ("shift", "None"): "cut_shift",
    }
    cutoff_style = cass_cutoffs[(job.sp.cutoff_style, job.sp.long_range_correction)]

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

    if any([abs(a.charge)>0.0 for a in structure.atoms]):
        charge_style = "ewald"
    else:
        charge_style = "none"

    Tmelt = 1000.0 * u.K

    T = job.sp.temperature * u.K
    P = job.sp.pressure * u.kPa

    species_list = [structure]
    cutoff = job.sp.r_cut * u.nm

    seedslist = [
        [7860904, 8601355],
        [5793508, 4173039],
        [4420642, 8720464],
        [8120272, 5850411],
        [7616664, 1492980],
        [6844679, 6087693],
        [2175335, 1317929],
        [9725500, 6331893],
        [4247127, 1385831],
        [2946981, 9870819],
        [8434295, 8017520],
        [8424221, 4595446],
        [8870203, 3009902],
        [4564019, 4788324],
        [3927152, 2536489],
        [3375750, 1798462],
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
                cutoff_style=cutoff_style,
                charge_style=charge_style,
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
                cutoff_style=cutoff_style,
                charge_style=charge_style,
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
                    cutoff_style=cutoff_style,
                    charge_style=charge_style,
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
                    cutoff_style=cutoff_style,
                    charge_style=charge_style,
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

        this_run = "equil_" + str(cycles_done)
        if not has_checkpoint(this_run):
            mc.run(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=equil_length,
                temperature=T,
                pressure=P,  # this line is ignored if ensemble isn't NPT
                properties=proplist,
                cutoff_style=cutoff_style,
                charge_style=charge_style,
                vdw_cutoff=cutoff,
                charge_cutoff=cutoff,
                run_name=this_run,
                prop_freq=prop_freq,
                coord_freq=coord_freq,
                units="sweeps",
                steps_per_sweep=N,
                seeds=seeds,
            )
        elif not check_complete(this_run):
            mc.restart(
                restart_from=get_last_checkpoint(this_run),
            )

        cycles_done += equil_length

        if ensemble == "GEMC-NVT":
            prpsuffix_liq = ".out.box1.prp"
            prpsuffix_vap = ".out.box2.prp"
            xyzsuffix_liq = ".out.box1.xyz"
            xyzsuffix_vap = ".out.box2.xyz"
            Hsuffix_liq = ".out.box1.H"
            Hsuffix_vap = ".out.box2.H"
            merge_restart_prp(this_run + prpsuffix_vap)
        else:
            prpsuffix_liq = ".out.prp"
            xyzsuffix_liq = ".out.xyz"
            Hsuffix_liq = ".out.H"
        merge_restart_prp(this_run + prpsuffix_liq)

        t, g, Neff_max = detectEquilibration(
            np.loadtxt(this_run + prpsuffix_liq, usecols=5)
        )

        if ensemble == "GEMC-NVT":
            tvap, gvap, Neff_max_vap = detectEquilibration(
                np.loadtxt(this_run + prpsuffix_vap, usecols=5)
            )
            t = max(t, tvap)
        prior_run = this_run

        for i in range(2):
            if t >= equil_length * 3 / (prop_freq * 4):
                this_run = "equil_" + str(cycles_done)
                if not has_checkpoint(this_run):
                    mc.restart(
                        restart_from=prior_run,
                        run_name=this_run,
                        run_type="equilibration",
                        total_run_length=cycles_done + equil_length,
                    )
                elif not check_complete(this_run):
                    mc.restart(
                        restart_from=get_last_checkpoint(this_run),
                    )
                cycles_done += equil_length
                merge_restart_prp(this_run + prpsuffix_liq)
                t, g, Neff_max = detectEquilibration(
                    np.loadtxt(this_run + prpsuffix_liq, usecols=5)
                )

                if ensemble == "GEMC-NVT":
                    merge_restart_prp(this_run + prpsuffix_vap)
                    tvap, gvap, Neff_max_vap = detectEquilibration(
                        np.loadtxt(this_run + prpsuffix_vap, usecols=5)
                    )
                    t = max(t, tvap)

                for rmtgt in list(Path(".").glob("equil_*.[xH]*")):
                    os.remove(rmtgt)
                prior_run = this_run

            else:
                break

        if has_checkpoint("prod"):
            mc.restart(
                restart_from=get_last_checkpoint("prod"),
            )
        else:
            mc.restart(
                restart_from=prior_run,
                run_name="prod",
                run_type="production",
                total_run_length=cycles_done + prod_length,
            )
        merge_restart_prp("prod" + prpsuffix_liq)
        merge_restart_traj("prod" + xyzsuffix_liq)
        merge_restart_traj("prod" + Hsuffix_liq)
        if ensemble == "GEMC-NVT":
            merge_restart_prp("prod" + prpsuffix_vap)
            merge_restart_traj("prod" + xyzsuffix_vap)
            merge_restart_traj("prod" + Hsuffix_vap)


@Project.operation
@Project.pre(cassandra_complete)
@Project.post(lambda job: "mean_energy_box1" in job.document)
def statistics(job):
    """Compute statistical quantities for each job."""
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

    if job.sp.ensemble == "GEMC-NVT":

        box1 = "prod.out.box1.prp"
        box2 = "prod.out.box2.prp"

        data_box1 = pd.read_table(
            job.fn(box1), skiprows=3, sep="\s+", names=proplist
        )
        data_box2 = pd.read_table(
            job.fn(box2), skiprows=3, sep="\s+", names=proplist
        )

        job.document.mean_energy_box1 = data_box1["energy_total"].mean()
        job.document.mean_energy_box2 = data_box2["energy_total"].mean()

    else:

        box1 = "prod.out.prp"

        data_box1 = pd.read_table(
            job.fn(box1), skiprows=3, sep="\s+", names=proplist
        )

        job.document.mean_energy_box1 = data_box1["energy_total"].mean()
        job.document.mean_density_box1 = (
            data_box1["mass_density"].mean() * 0.001
        )


def prp2txt(prp_path, txt_path, T):
    """Convert Cassandra prp file to project-standard txt file."""
    import numpy as np

    prp_array = np.loadtxt(prp_path, usecols=(1, 4, 6))
    potential_energy = prp_array[:, 0]
    pressure = prp_array[:, 1] * 100.0
    density = prp_array[:, 2] * 0.001
    n_points = len(density)
    kinetic_energy = np.zeros(n_points, dtype=int)
    temperature = np.full(n_points, T)
    timestep = np.arange(10, (n_points + 1) * 10, 10)
    prpdf = pd.DataFrame(
        {
            "timestep": timestep,
            "potential_energy": potential_energy,
            "kinetic_energy": kinetic_energy,
            "pressure": pressure,
            "temperature": temperature,
            "density": density,
        }
    )
    prpdf.to_csv(txt_path, sep=" ", index=False)


@Project.label
def output_processed(job):
    """Check whether Cassandra's output has been processed."""
    return os.path.exists(job.fn("log-npt.txt")) or os.path.exists(
        job.fn("log-vapor.txt")
    )


@Project.operation
@Project.pre(cassandra_complete)
@Project.post(output_processed)
def process_output(job):
    """Convert Cassandra trajectories to gsd format and convert Cassandra property output to project-standard txt files."""
    import foyer
    import mbuild as mb
    import numpy as np

    from reproducibility_project.src.molecules.system_builder import (
        get_molecule,
    )
    from reproducibility_project.src.utils.forcefields import load_ff
    from reproducibility_project.src.utils.trajectory_conversion import (
        cassandra2gsd,
    )

    species_list = [load_ff(job.sp.forcefield_name).apply(get_molecule(job.sp))]
    if job.sp.ensemble == "GEMC-NVT":
        liqbox = "prod.out.box1"
        vapbox = "prod.out.box2"
        liq_H = job.fn(liqbox + ".H")
        liq_xyz = job.fn(liqbox + ".xyz")
        vap_H = job.fn(vapbox + ".H")
        vap_xyz = job.fn(vapbox + ".xyz")
        liq_prp = job.fn(liqbox + ".prp")
        vap_prp = job.fn(vapbox + ".prp")
        f_list = [liq_H, liq_xyz, vap_H, vap_xyz, liq_prp, vap_prp]
        if not all([os.path.exists(f) for f in f_list]):
            return
        cassandra2gsd(
            h_path=liq_H,
            xyz_path=liq_xyz,
            gsd_path=job.fn("trajectory-liquid.gsd"),
            species_list=species_list,
        )
        cassandra2gsd(
            h_path=vap_H,
            xyz_path=vap_xyz,
            gsd_path=job.fn("trajectory-vapor.gsd"),
            species_list=species_list,
        )
        prp2txt(liq_prp, job.fn("log-liquid.txt"), job.sp.temperature)
        prp2txt(vap_prp, job.fn("log-vapor.txt"), job.sp.temperature)
    else:
        nptbox = "prod.out"
        npt_H = job.fn(nptbox + ".H")
        npt_xyz = job.fn(nptbox + ".xyz")
        npt_prp = job.fn(nptbox + ".prp")
        f_list = [npt_H, npt_xyz, npt_prp]
        if not all([os.path.exists(f) for f in f_list]):
            return
        cassandra2gsd(
            h_path=npt_H,
            xyz_path=npt_xyz,
            gsd_path=job.fn("trajectory-npt.gsd"),
            species_list=species_list,
        )
        prp2txt(npt_prp, job.fn("log-npt.txt"), job.sp.temperature)


if __name__ == "__main__":
    pr = Project()
    pr.main()

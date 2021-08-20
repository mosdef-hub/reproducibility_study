"""Setup for signac, signac-flow, signac-dashboard for this study."""
# import foyer
import os
import pathlib

import flow
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


@Project.label
def is_cassandra(job):
    return job.sp.simulation_engine == "cassandra"

@Project.operation
@Project.post(cassandra_complete)
@directives(omp_num_threads=4)
def run_cassandra(job):
    ""

    import mbuild as mb
    import mosdef_cassandra as mc
    import foyer
    import unyt as u
    from mbuild.formats.xyz import read_xyz

    from reproducibility_project.src.molecules.system_builder import construct_system
    from reproducibility_project.src.molecules.methane_ua import MethaneUA
    from reproducibility_project.src.molecules.pentane_ua import PentaneUA


    compound_dict = {
            "methaneUA": MethaneUA(),
            "pentaneUA": PentaneUA(),
            "benzeneUA": BenzeneUA(),
            "waterSPC/E": mb.lib.molecules.water.WaterSPC(),
            "ethanolAA": mb.load("CCO", smiles=True),
            }
    molecule = job.sp.molecule

    compound = compound_dict[molecule]

    ensemble = job.sp.production_ensemble

    Nliq = job.sp.N_liquid
    Nvap = job.sp.N_vap
    N = Nliq+Nvap

    filled_boxes = construct_system(job.sp)

    liqbox_filled = filled_boxes[0]
    if ensemble == "GEMC-NVT":
        vapbox_filled = filled_boxes[1]

    ffname = job.sp.forcefield
    ff = get_ff(ffname) # make this function exist
    structure = ff.apply(compound)

    T = job.sp.temperature * u.K
    P = job.sp.pressure * u.kPa


    species_list = [structure]
    cutoff_dict = {
            "Trappe_UA": 10.0,
            "spce": 9.0,
            "oplsaa": 14.0
    cutoff = cutoff_dict[ffname] * u.Angstrom

    seedslist = [
        [576243192, 277516412],
        [276640023, 802773011],
        [911197026, 307764101],
        [774491056, 422678105],
        [723488469, 111294037],
        ]

    seeds = seedslist[job.sp.replica]


    cbmc_n_ins = 12
    cbmc_n_dihed = 50

    prop_freq = 10
    coord_freq = 10
    if ff_name == "oplsaa":
        comb_rule = 'geometric'
    else:
        comb_rule = 'lb'


    proplist = [
            "energy_total",
            "volume",
            "nmols",
            "pressure"
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
            "enthalpy"]




    meltsystem_liq = mc.System([liqbox_filled], species_list)

    nvtmoves = mc.MoveSet('nvt', species_list)
    nvtmoves.prob_rotate = p_rotate
    nvtmoves.prob_translate = p_translate
    nvtmoves.prob_regrow = p_regrow

    moveset = mc.MoveSet(ensemble, species_list)
    moveset.prob_volume = p_volume
    moveset.prob_translate = p_translate
    moveset.prob_rotate = p_rotate
    moveset.prob_regrow = p_regrow
    moveset.cbmc_n_ins = cbmc_n_ins
    moveset.cbmc_n_dihed = cbmc_n_dihed





    mc.run(
            system=meltsystem_liq,
            moveset=nvtmoves,
            run_type="equilibration",
            run_length=5000,
            temperature=1000,
            pressure=P,
            properties=proplist,
            cutoff_style="cut",
            vdw_cutoff=cutoff,
            charge_cutoff=cutoff,
            run_name="nvt_melt",
            prop_freq=prop_freq,
            coord_freq=5000,
            units='sweeps',
            steps_per_sweep=Nliq,
            seeds=seeds
            )

    nvtendbox_liq = read_xyz("nvt_melt.out.xyz")

    nvtendbox_liq.box = liqbox_filled.box

    nvtsystem_liq = mc.System([nvtendbox_liq], species_list)

    mc.run(
            system=nvtsystem_liq,
            moveset=nvtmoves,
            run_type="equilibration",
            run_length=5000,
            temperature=T,
            pressure=P,
            properties=proplist,
            cutoff_style="cut",
            vdw_cutoff=cutoff,
            charge_cutoff=cutoff,
            run_name="nvt_equil",
            prop_freq=prop_freq,
            coord_freq=5000,
            units='sweeps',
            steps_per_sweep=Nliq,
            seeds=seeds
            )

    nvtendbox_liq = read_xyz("nvt_equil.out.xyz")

    nvtendbox_liq.box = liqbox_filled.box
    boxlist = [nvtendbox_liq]

    if ensemble == "GEMC-NVT":
        meltsystem_vap = mc.System([vapbox_filled], species_list)
        mc.run(
                system=meltsystem_vap,
                moveset=nvtmoves,
                run_type="equilibration",
                run_length=5000,
                temperature=1000,
                pressure=P,
                properties=proplist,
                cutoff_style="cut",
                vdw_cutoff=cutoff,
                charge_cutoff=cutoff,
                run_name="nvt_melt_vap",
                prop_freq=prop_freq,
                coord_freq=5000,
                units='sweeps',
                steps_per_sweep=Nvap,
                seeds=seeds
                )
        meltendbox_vap = read_xyz("nvt_melt_vap.out.xyz")

        meltendbox_vap.box = vapbox_filled.box
        nvtsystem_vap = mc.System([meltendbox_vap], species_list)
        mc.run(
                system=nvtsystem_vap,
                moveset=nvtmoves,
                run_type="equilibration",
                run_length=5000,
                temperature=T,
                pressure=P,
                properties=proplist,
                cutoff_style="cut",
                vdw_cutoff=cutoff,
                charge_cutoff=cutoff,
                run_name="nvt_equil_vap",
                prop_freq=prop_freq,
                coord_freq=5000,
                units='sweeps',
                steps_per_sweep=Nvap,
                seeds=seeds
                )
        nvtendbox_vap = read_xyz("nvt_equil_vap.out.xyz")

        nvtendbox_vap.box = vapbox_filled.box
        boxlist.append(nvtendbox_vap)
        moveset.prob_swap = p_swap


    system = mc.System(boxlist, species_list)

    mc.run(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=50000,
            temperature=T,
            pressure=P,
            properties=proplist,
            cutoff_style="cut",
            vdw_cutoff=cutoff,
            charge_cutoff=cutoff,
            run_name="equil",
            prop_freq=prop_freq,
            coord_freq=coord_freq,
            units='sweeps',
            steps_per_sweep=N,
            seeds=seeds
            )

    mc.restart(
            restart_from="equil",
            run_name="prod",
            run_type="production",
            total_run_length=170000
            )

if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()

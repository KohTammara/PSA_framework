/* 
 * 'PSA_framework' - A Nextflow pipeline for the analysis of protein structure either designed or mutated
 * 
 * Tammara Poa Siang Koh 
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2


/*
 * Define the parameters
 */ 


 log.info """\
C A L L I N G S  -  PSA    v 1 
================================
list of structures (threader/FBB)       : $params.list_of_structs
residue file (FBB)                      : $params.resf
sequence file (threader)                : $params.sequence

mutations (threader/STRUM)              : $params.mutation
name for threading mover                : $params.name
Pack neighbor (threader)                : $params.pack_neighbors
Neightbor distane (threader)            : $params.neighbor_dis
score function name (threader)          : $params.scorefxn
weights (threader)                      : $params.weights
skip unknown mutants (threader)         : $params.skip_unknown_mutant
pack round (threader)                   : $params.pack_round
sequence mode (threader)                : $params.sequence_mode
xml template (threader)                 : $params.template

"""
/* 
 * Import modules and alias for reuse
 */
include { 
    ROSETTA_FIXBB;
    ROSETTA_THREADER;
    MAESTRO_XML;
    MAESTRO;
    CUTANDMUTATE;
    CREATEXML;
    STANDARD_MD as GROMACS_MT_FBB;
    STANDARD_MD as GROMACS_MT_THREADER;
    GROMACS_MT_PMX;
    STANDARD_MD as GROMACS_WT;
    PMX_PREP_MUTANT;
    GRO_PREP_MUTANT as GRO_PREP_MUTANT_FBB;
    GRO_PREP_MUTANT as GRO_PREP_MUTANT_THR;
    GRO_EQUILIBRIUM as GRO_EQUILIBRIUM_FOR;
    GRO_EQUILIBRIUM as GRO_EQUILIBRIUM_FOR_FBB;
    GRO_EQUILIBRIUM as GRO_EQUILIBRIUM_FOR_THR;
    GRO_EQUILIBRIUM as GRO_EQUILIBRIUM_REV;
    GRO_EQUILIBRIUM as GRO_EQUILIBRIUM_REV_FBB;
    GRO_EQUILIBRIUM as GRO_EQUILIBRIUM_REV_THR;
    GRO_NON_EQUILIBRIUM as GRO_NON_EQUILIBRIUM_FOR;
    GRO_NON_EQUILIBRIUM as GRO_NON_EQUILIBRIUM_FOR_FBB;
    GRO_NON_EQUILIBRIUM as GRO_NON_EQUILIBRIUM_FOR_THR;
    GRO_NON_EQUILIBRIUM_UNFOLDED as GRO_NON_EQUILIBRIUM_FOR_UNF;
    GRO_NON_EQUILIBRIUM as GRO_NON_EQUILIBRIUM_REV;
    GRO_NON_EQUILIBRIUM as GRO_NON_EQUILIBRIUM_REV_FBB;
    GRO_NON_EQUILIBRIUM as GRO_NON_EQUILIBRIUM_REV_THR;
    GRO_NON_EQUILIBRIUM_UNFOLDED as GRO_NON_EQUILIBRIUM_REV_UNF;
    FREE_ENERGY_EST as FEE_FBB;
    FREE_ENERGY_EST as FEE_THR;
    FREE_ENERGY_EST as FEE_PMX;
    FREE_ENERGY_EST as FEE_PMX_UNF;
    CREATE_TRIPEPTIDE;
    READ_TRIPEPTIDE_FILES;
    GRO_PREP_TRIPEPTIDE;
    GRO_EQUILIBRIUM_UNFOLDED as GRO_EQUILIBRIUM_UNFOLDED_FOR;
    GRO_EQUILIBRIUM_UNFOLDED as GRO_EQUILIBRIUM_UNFOLDED_REV;
} from './modules.nf' 

workflow.onComplete {
    println ""
    println "Workflow info:"
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

//workflow to execute rosetta fbb
workflow rosy_fbb {
    take:
    pdb
    resf

    main:
    ROSETTA_FIXBB(pdb, resf)

    emit:
    pdb = ROSETTA_FIXBB.out[0]
    score = ROSETTA_FIXBB.out[1]
    log_file = ROSETTA_FIXBB.out[2]
    mutation_name = ROSETTA_FIXBB.out[3]
}

//workflow to create input for threading with rosetta
workflow rosy_threader_input {
    take:
    sequence
    mutation

    main:
    CUTANDMUTATE(sequence, mutation)

    emit:
    cut_seq = CUTANDMUTATE.out[0]
    start_position = CUTANDMUTATE.out[1]
}

//workflow to execute rosetta threading
workflow rosy_threader {
    take:
    pdb
    cut_sequence
    start_position
    mutation

    main:
    xml = CREATEXML(params.name, cut_sequence, params.pack_round, params.sequence_mode, params.skip_unknown_mutant, params.scorefxn, start_position, params.neighbor_dis, params.pack_neighbors, params.weights, params.template)
    ROSETTA_THREADER(pdb, xml, mutation)

    emit:
    pdb = ROSETTA_THREADER.out[0]
    score = ROSETTA_THREADER.out[1]
    mutation_name = ROSETTA_THREADER.out[2]
}

//workflow to execute maestro
workflow maestro {
    take:
    mut
    xml

    main:
    MAESTRO(params.effiles, params.council, mut, xml)

    emit:
    results_csv = MAESTRO.out[0]
}

//workflow for free energy forward transition of the pmx mutated folded state
workflow pmx_free_energy_forward {
    take:
    ions_pdb
    newtop
    posre_itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_FOR(ions_pdb, newtop, posre_itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_FOR.output[0]
    equi_tpr = GRO_EQUILIBRIUM_FOR.output[1]
    name = GRO_EQUILIBRIUM_FOR.output[4]
    GRO_NON_EQUILIBRIUM_FOR(equi_trr,equi_tpr,newtop, nonequil, name, "for")

    emit:
    forward = GRO_NON_EQUILIBRIUM_FOR.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_FOR.out[1]

}

//workflow for free energy forward transition of the pmx unfolded state (tripeptide)
workflow pmx_free_energy_forward_unfolded {
    take:
    ions_pdb
    newtop
    posre_itp
    itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_UNFOLDED_FOR(ions_pdb, newtop, posre_itp, itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_UNFOLDED_FOR.output[0]
    equi_tpr = GRO_EQUILIBRIUM_UNFOLDED_FOR.output[1]
    name = GRO_EQUILIBRIUM_UNFOLDED_FOR.output[4]
    top = GRO_EQUILIBRIUM_UNFOLDED_FOR.output[5]
    n_posre_itp = GRO_EQUILIBRIUM_UNFOLDED_FOR.output[6]
    n_itp = GRO_EQUILIBRIUM_UNFOLDED_FOR.output[7]
    GRO_NON_EQUILIBRIUM_FOR_UNF(equi_trr, equi_tpr, top, n_posre_itp, n_itp, nonequil, name, "for")

    emit:
    forward = GRO_NON_EQUILIBRIUM_FOR_UNF.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_FOR_UNF.out[1]

}

//workflow for free energy reverse transition of the pmx mutated folded state
workflow pmx_free_energy_reverse {
    take:
    ions_pdb
    newtop
    posre_itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_REV(ions_pdb, newtop, posre_itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_REV.output[0]
    equi_tpr = GRO_EQUILIBRIUM_REV.output[1]
    GRO_NON_EQUILIBRIUM_REV(equi_trr,equi_tpr,newtop, nonequil, name, "rev")

    emit:
    reverse =  GRO_NON_EQUILIBRIUM_REV.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_REV.out[1]
}

//workflow for free energy reverse transition of the pmx unfolded state (tripeptide)
workflow pmx_free_energy_reverse_unfolded {
    take:
    ions_pdb
    newtop
    posre_itp
    itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_UNFOLDED_REV(ions_pdb, newtop, posre_itp, itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_UNFOLDED_REV.output[0]
    equi_tpr = GRO_EQUILIBRIUM_UNFOLDED_REV.output[1]
    name = GRO_EQUILIBRIUM_UNFOLDED_REV.output[4]
    top = GRO_EQUILIBRIUM_UNFOLDED_REV.output[5]
    n_posre_itp = GRO_EQUILIBRIUM_UNFOLDED_REV.output[6]
    n_itp = GRO_EQUILIBRIUM_UNFOLDED_REV.output[7]
    GRO_NON_EQUILIBRIUM_REV_UNF(equi_trr, equi_tpr, top, n_posre_itp, n_itp, nonequil, name, "rev")

    emit:
    reverse =  GRO_NON_EQUILIBRIUM_REV_UNF.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_REV_UNF.out[1]
}

//workflow that executes the whole free energy estimation process for the pmx folded state
workflow free_energy_folded {
    take:
    ions_pdb
    newtop
    posre_itp
    f_enmin
    f_equil
    f_npt
    f_nonequil
    name
    r_enmin
    r_equil
    r_npt
    r_nonequil

    main:
    pmx_free_energy_forward(ions_pdb, newtop, posre_itp, f_enmin, f_equil, f_npt, f_nonequil, name)
    pmx_free_energy_reverse(ions_pdb, newtop, posre_itp, r_enmin, r_equil, r_npt, r_nonequil, name)
    FEE_PMX(pmx_free_energy_forward.out.forward, pmx_free_energy_reverse.out.reverse, "folded", pmx_free_energy_forward.out.name)

    emit:
    FEE_PMX.out[0]
}

//workflow that executes the whole free energy estimation process for the pmx unfolded state
workflow free_energy_unfolded {
    take:
    ions_pdb
    newtop
    posre_itp
    itp
    f_enmin
    f_equil
    f_npt
    f_nonequil
    name
    r_enmin
    r_equil
    r_npt
    r_nonequil

    main:
    pmx_free_energy_forward_unfolded(ions_pdb, newtop, posre_itp, itp, f_enmin, f_equil, f_npt, f_nonequil, name)
    pmx_free_energy_reverse_unfolded(ions_pdb, newtop, posre_itp, itp, r_enmin, r_equil, r_npt, r_nonequil, name)
    FEE_PMX_UNF(pmx_free_energy_forward_unfolded.out.forward, pmx_free_energy_reverse_unfolded.out.reverse, "unfolded", pmx_free_energy_forward_unfolded.out.name)

    emit:
    FEE_PMX_UNF.out[0]
}

//workflow for free energy forward transition of the rosetta fbb mutated folded state
workflow free_energy_forward_fbb {
    take:
    ions_pdb
    newtop
    posre_itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_FOR_FBB(ions_pdb, newtop, posre_itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_FOR_FBB.output[0]
    equi_tpr = GRO_EQUILIBRIUM_FOR_FBB.output[1]
    name = GRO_EQUILIBRIUM_FOR_FBB.output[4]
    GRO_NON_EQUILIBRIUM_FOR_FBB(equi_trr,equi_tpr,newtop, nonequil, name, "for")

    emit:
    forward = GRO_NON_EQUILIBRIUM_FOR_FBB.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_FOR_FBB.out[1]

}

//workflow for free energy reverse transition of the rosetta fbb mutated folded state
workflow free_energy_reverse_fbb {
    take:
    ions_pdb
    newtop
    posre_itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_REV_FBB(ions_pdb, newtop, posre_itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_REV_FBB.output[0]
    equi_tpr = GRO_EQUILIBRIUM_REV_FBB.output[1]
    GRO_NON_EQUILIBRIUM_REV_FBB(equi_trr,equi_tpr,newtop, nonequil, name, "rev")

    emit:
    reverse =  GRO_NON_EQUILIBRIUM_REV_FBB.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_REV_FBB.out[1]
}

//workflow that executes the whole free energy estimation process for the rosetta fbb folded state
workflow free_energy_folded_fbb {
    take:
    ions_pdb
    newtop
    posre_itp
    f_enmin
    f_equil
    f_npt
    f_nonequil
    name
    r_enmin
    r_equil
    r_npt
    r_nonequil

    main:
    free_energy_forward_fbb(ions_pdb, newtop, posre_itp, f_enmin, f_equil, f_npt, f_nonequil, name)
    free_energy_reverse_fbb(ions_pdb, newtop, posre_itp, r_enmin, r_equil, r_npt, r_nonequil, name)
    FEE_FBB(free_energy_forward_fbb.out.forward, free_energy_reverse_fbb.out.reverse, "folded", free_energy_forward_fbb.out.name)

    emit:
    FEE_FBB.out[0]
}

//workflow for free energy forward transition of the rosetta threaded mutated folded state
workflow free_energy_forward_thr {
    take:
    ions_pdb
    newtop
    posre_itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_FOR_THR(ions_pdb, newtop, posre_itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_FOR_THR.output[0]
    equi_tpr = GRO_EQUILIBRIUM_FOR_THR.output[1]
    name = GRO_EQUILIBRIUM_FOR_THR.output[4]
    GRO_NON_EQUILIBRIUM_FOR_THR(equi_trr,equi_tpr,newtop, nonequil, name, "for")

    emit:
    forward = GRO_NON_EQUILIBRIUM_FOR_THR.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_FOR_THR.out[1]

}

//workflow for free energy reverse transition of the rosetta threaded mutated folded state
workflow free_energy_reverse_thr {
    take:
    ions_pdb
    newtop
    posre_itp
    enmin
    equil
    npt
    nonequil
    name

    main:
    GRO_EQUILIBRIUM_REV_THR(ions_pdb, newtop, posre_itp, enmin, equil, npt, name)
    equi_trr = GRO_EQUILIBRIUM_REV_THR.output[0]
    equi_tpr = GRO_EQUILIBRIUM_REV_THR.output[1]
    GRO_NON_EQUILIBRIUM_REV_THR(equi_trr,equi_tpr,newtop, nonequil, name, "rev")

    emit:
    reverse =  GRO_NON_EQUILIBRIUM_REV_THR.out[0].collect()
    name = GRO_NON_EQUILIBRIUM_REV_THR.out[1]
}

//workflow that executes the whole free energy estimation process for the rosetta threaded folded state
workflow free_energy_folded_thr {
    take:
    ions_pdb
    newtop
    posre_itp
    f_enmin
    f_equil
    f_npt
    f_nonequil
    name
    r_enmin
    r_equil
    r_npt
    r_nonequil

    main:
    free_energy_forward_thr(ions_pdb, newtop, posre_itp, f_enmin, f_equil, f_npt, f_nonequil, name)
    free_energy_reverse_thr(ions_pdb, newtop, posre_itp, r_enmin, r_equil, r_npt, r_nonequil, name)
    FEE_THR(free_energy_forward_thr.out.forward, free_energy_reverse_thr.out.reverse, "folded", free_energy_forward_thr.out.name)

    emit:
    FEE_THR.out[0]
}


workflow {
    if (params.pmx == true) {
        //read tsv mutation data for pmx folded 
        pmx_mutants_folded = Channel.fromPath(params.mutation_file_folded).splitCsv(header: true, sep: '\t').map { row -> tuple(row.index, row.mutant) }
        PMX_PREP_MUTANT(params.pdb, pmx_mutants_folded, params.ff_name, params.genion_mdp)
        newtop = PMX_PREP_MUTANT.output[1]
        ions_pdb = PMX_PREP_MUTANT.output[6]
        posre_itp = PMX_PREP_MUTANT.output[7]
        mutation_name = PMX_PREP_MUTANT.output[8]
        free_energy_folded(ions_pdb, newtop, posre_itp, params.f_enmin_mdp, params.f_equil_mdp, params.f_npt_mdp, params.f_nonequil_mdp, mutation_name, params.r_enmin_mdp, params.r_equil_mdp, params.r_npt_mdp, params.r_nonequil_mdp)
        
    }
    if (params.existing_tripeptide_files == true) {
        tripep_dir = Channel
            .fromPath("${params.tripeptide_files}/tripep_*", type: 'dir')
            .ifEmpty { error "No tripeptide folders found in ${params.tripeptide_files}" }
            .map { folder -> 
                folder 
            }
        files = READ_TRIPEPTIDE_FILES(tripep_dir)
        prep_files = GRO_PREP_TRIPEPTIDE(files,params.genion_mdp,params.ff_name)
        free_energy_unfolded(prep_files[5], prep_files[0], prep_files[6], prep_files[7], params.f_enmin_mdp, params.f_equil_mdp, params.f_npt_mdp, params.f_nonequil_mdp, prep_files[8], params.r_enmin_mdp, params.r_equil_mdp, params.r_npt_mdp, params.r_nonequil_mdp)

    }

    if (params.rosetta_fbb == true) {
        resf_folded = Channel.fromPath(params.resf_folded)
        rosy_fbb(params.list_of_structs, resf_folded)
        //free energy folded
        GRO_PREP_MUTANT_FBB(rosy_fbb.out.pdb, rosy_fbb.out.mutation_name, params.ff_name, params.genion_mdp)
        topol = GRO_PREP_MUTANT_FBB.output[0]
        ions_pdb = GRO_PREP_MUTANT_FBB.output[5]
        posre_itp = GRO_PREP_MUTANT_FBB.output[6]
        mutation_name = GRO_PREP_MUTANT_FBB.output[7]
        free_energy_folded_fbb(ions_pdb, topol, posre_itp, params.f_enmin_mdp, params.f_equil_mdp, params.f_npt_mdp, params.f_nonequil_mdp,mutation_name, params.r_enmin_mdp, params.r_equil_mdp, params.r_npt_mdp, params.r_nonequil_mdp)
        
    }

    if (params.rosetta_threader == true) {
        mutation = Channel.fromPath(params.mutation_info)
        rosy_threader_input(params.sequence, mutation)
        rosy_threader(params.list_of_structs,rosy_threader_input.out.cut_seq,rosy_threader_input.out.start_position, mutation)

        GRO_PREP_MUTANT_THR(rosy_threader.out.pdb, rosy_threader.out.mutation_name, params.ff_name, params.genion_mdp)
        topol = GRO_PREP_MUTANT_THR.output[0]
        ions_pdb = GRO_PREP_MUTANT_THR.output[5]
        posre_itp = GRO_PREP_MUTANT_THR.output[6]
        mutation_name = GRO_PREP_MUTANT_THR.output[7]
        free_energy_folded_thr(ions_pdb, topol, posre_itp, params.f_enmin_mdp, params.f_equil_mdp, params.f_npt_mdp, params.f_nonequil_mdp,mutation_name,params.r_enmin_mdp, params.r_equil_mdp, params.r_npt_mdp, params.r_nonequil_mdp)
    }

    if (params.maestro == true) {
        xml = MAESTRO_XML(params.effiles, params.council, params.path_to_pdb, params.prefix_maestro_out, params.postfix_maestro_out, params.to_lower_maestro_out, params.to_upper_maestro_out, params.bu_maestro) 
        maestro(params.mutationfile, xml)
    }

    if (params.gromacs == true) {
        if (params.rosetta_fbb == true) {
            GROMACS_MT_FBB(rosy_fbb.out.pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9, params.G10)
            GROMACS_WT(params.pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9, params.G10)
        }
        if (params.rosetta_threader == true) {
            GROMACS_MT_THREADER(rosy_threader.out.pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9, params.G10)
            GROMACS_WT(params.pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9, params.G10)
        }
        if (params.rosetta_threader == false && params.rosetta_threader == false) {
            GROMACS_WT(params.pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9, params.G10)
        }
        if (params.pmx == true) {
            GROMACS_MT_PMX(PMX_PREP_MUTANT.output[6],PMX_PREP_MUTANT.output[1],PMX_PREP_MUTANT.output[7], params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9, params.G10)
        }
    }

    // publishDir '/results', mode: 'copy', overwrite: false
    // Rosetta fbb execution
    // pdb = Channel.fromPath(params.list_of_structs)
    // resf = Channel.fromPath(params.resf)
    // rosy_fbb(pdb, resf)
    // // rosy_fbb.out.pdb.view()
    // // rosy_fbb.out.score.view()
    // // rosy_fbb.out.log_file.view()
    

    //Rosetta threader execution
    // sequence = Channel.fromPath(params.sequence)
    // mutation = Channel.fromPath(params.mutation_info)
    // // sequence.view()
    // // mutation.view()
    // rosy_threader_input(sequence, mutation)
    // // rosy_threader_input.out.cut_seq.view()
    // rosy_threader(pdb,rosy_threader_input.out.cut_seq,rosy_threader_input.out.start_position)

    //MAESTRO execution
    // maestro(pdb)
    // // maestro.out.results_csv.view()


    // GROMACS execution FBB
    // GROMACS_MT_FBB(rosy_fbb_pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
    // GROMACS_MT_THREADER(rosy_threader_pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
    // GROMACS_WT(pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
}

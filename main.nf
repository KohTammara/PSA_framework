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
 * Import modules 
 */
include { 
    ROSETTA_FIXBB;
    ROSETTA_THREADER;
    MAESTRO_XML;
    MAESTRO;
    CUTANDMUTATE;
    CREATEXML;
    GROMACS_MT_FBB;
    GROMACS_MT_THREADER;
    GROMACS_WT;
    PMX_PREP_MUTANT;
    GRO_PREP_MUTANT;
    GRO_EQUILIBRIUM;
    // GRO_NON_EQUILIBRIUM;
} from './modules.nf' 

workflow.onComplete {
    println ""
    println "Workflow info:"
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

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
}

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

workflow rosy_threader {
    take:
    pdb
    cut_sequence
    start_position
    mutation

    main:
    xml = CREATEXML(params.name, cut_sequence, params.sequence_mode, params.pack_round, params.skip_unknown_mutant, params.scorefxn, start_position, params.neighbor_dis, params.pack_neighbors, params.weights, params.template)
    ROSETTA_THREADER(pdb, xml, mutation)

    emit:
    pdb = ROSETTA_THREADER.out[0]
    score = ROSETTA_THREADER.out[1]
}

workflow maestro {
    take:
    mut
    xml

    main:
    MAESTRO(params.effiles, params.council, mut, xml)

    emit:
    results_csv = MAESTRO.out[0]
}


workflow {
    if (params.pmx == true) {
        PMX_PREP_MUTANT(params.pdb, params.res_number, params.mutant_res, params.pmx_forcefield)
        newtop = PMX_PREP_MUTANT.output[1]
        ions_pdb = PMX_PREP_MUTANT.output[6]
        posre_itp = PMX_PREP_MUTANT.output[7]
        GRO_EQUILIBRIUM(ions_pdb, newtop, posre_itp)
        equi_trr = GRO_EQUILIBRIUM.output[0]
        equi_tpr = GRO_EQUILIBRIUM.output[1]
        newtop = GRO_EQUILIBRIUM.outPUT[4]
        GRO_NON_EQUILIBRIUM(equi_trr,equi_tpr,newtop)
    }

    if (params.rosetta_fbb == true) {
        // pdb = Channel.fromPath(params.list_of_structs)
        resf = Channel.fromPath(params.resf)
        resf.view()
        rosy_fbb(params.list_of_structs, resf)
    }

    if (params.rosetta_threader == true) {
        // pdb = Channel.fromPath(params.list_of_structs)
        // sequence = Channel.fromPath(params.sequence)
        mutation = Channel.fromPath(params.mutation_info)
        rosy_threader_input(params.sequence, mutation)
        rosy_threader(params.list_of_structs,rosy_threader_input.out.cut_seq,rosy_threader_input.out.start_position, mutation)
    }

    if (params.maestro == true) {
        xml = MAESTRO_XML(params.effiles, params.council, params.path_to_pdb, params.prefix_maestro_out, params.postfix_maestro_out, params.to_lower_maestro_out, params.bu_maestro) 
        maestro(params.mutationfile, xml)
    }

    if (params.gromacs == true) {
        if (params.rosetta_fbb == true) {
            GROMACS_MT_FBB(rosy_fbb.out.pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
            GROMACS_WT(pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
        }
        if (params.rosetta_threader == true) {
            GROMACS_MT_THREADER(rosy_threader.out.pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
            GROMACS_WT(pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
        }
        if ((params.rosetta_threader == false && params.rosetta_threader == false)||(params.rosetta_threader == true || params.rosetta_threader == true)) {
            GROMACS_WT(pdb, params.ions_mdp, params.em_mdp, params.nvt_mdp, params.npt_mdp, params.md_mdp, params.G1, params.G2, params.G3, params.G4, params.G5, params.G6,params.G7, params.G8, params.G9)
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

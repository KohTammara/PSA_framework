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
 * Define the default parameters
 */ 


 log.info """\
C A L L I N G S  -  PSA    v 1 
================================
list of structures (threader/FBB)       : $params.list_of_structs
residue file (FBB)                      : $params.resf
sequence file (threader)                : $params.sequence
start position of sequence (threader)   : $params.start_position
end position of sequence (threader)     : $params.end_position
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
    MAESTRO;
    CUTANDMUTATE;
    CREATEXML;
} from './modules.nf' 

/* 
 * main pipeline logic
 */
// workflow fbb {
//     //Rosetta fbb execution
//     pdb = Channel.fromPath( params.list_of_structs )
//     ROSETTA_FIXBB(pdb, params.resf)
// }
// workflow thread {
//     //Rosetta threader execution
//     seq = CUTANDMUTATE(params.sequence, params.start_position, params.end_position, params.mutation)
//     xml = CREATEXML(params.name, seq, params.sequence_mode, params.pack_round, params.skip_unknown_mutant, params.scorefxn, params.start_position, params.neighbor_dis, params.pack_neighbors, params.weights, params.template)
//     ROSETTA_THREADER(pdb, xml)
// }
workflow.onComplete {
    println ""
    println "Workflow info:"
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow {
    //Rosetta fbb execution
    pdb = Channel.fromPath( params.list_of_structs )
    ROSETTA_FIXBB(pdb, params.resf)

    //Rosetta threader execution
    seq = CUTANDMUTATE(params.sequence, params.start_position, params.end_position, params.mutation)
    xml = CREATEXML(params.name, seq, params.sequence_mode, params.pack_round, params.skip_unknown_mutant, params.scorefxn, params.start_position, params.neighbor_dis, params.pack_neighbors, params.weights, params.template)
    ROSETTA_THREADER(pdb, xml)

    // ddg monomer Rosetta [NO LONGER IN USE]
    // ROSETTA_DDG_PREMINIMIZATION(pdb)
    // ROSETTA_DDG(ROSETTA_DDG_PREMINIMIZATION(pdb), params.resf)

    //MAESTRO execution
    MAESTRO(params.effiles, params.council, pdb, params.mutation, params.chain, params.maestro_xml)
}

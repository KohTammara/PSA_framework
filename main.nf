/* 
 * 'PSA_framework' - A Nextflow pipeline for the analysis of protein structure either designed or mutated
 * 
 * Tammara Poa Siang Koh 
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

// // Ensure params.file is defined (adjust the path accordingly)
// params.file = file("params.config")

// // Load parameters from the config file
// params.file.text.readLines().each { line ->
//     def (key, value) = line.tokenize('=').collect { it.trim() }
//     println "Param $key: $value"
// }

/*
 * Define the default parameters
 */ 
 params.results    = "results"

 log.info """\
C A L L I N G S  -  PSA    v 1 
================================
results  : $params.results
"""
/* 
 * Import modules 
 */
include { 
    ROSETTA_FIXBB;
    ROSETTA_THREADER;
    CUTANDMUTATE;
    CREATEXML;
    SPLITPDB 
} from './modules.nf' 

/* 
 * main pipeline logic
 */
workflow {
    // pdb = SPLITPDB(params.list_of_structs)
    // pdb.view()
    // colllection = Channel.fromPath(params.list_of_structs).flatten()
    // ROSETTA_FIXBB(pdb, params.resf)


    //Rosetta fbb execution
    // params.list_pdb = "${workDir}/test/*.pdb"

    pdb = Channel.fromPath( params.list_of_structs )
    ROSETTA_FIXBB(pdb, params.resf)

    //Rosetta threader execution
    seq = CUTANDMUTATE(params.sequence, params.start_position, params.end_position, params.mutation)
    xml = CREATEXML(params.name, seq, params.sequence_mode, params.pack_round, params.skip_unknown_mutant, params.scorefxn, params.start_position, params.neighbor_dis, params.pack_neighbors, params.weights, params.template)
    ROSETTA_THREADER(pdb, xml)
}

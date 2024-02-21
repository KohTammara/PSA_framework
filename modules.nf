

// process STRUM {
// 	/*
// 	Adapt to read from file for multiple instances
// 	*/
// 	input:
// 	val mutation
// 	val input_type
// 	path path_to_in
	
// 	output:
//     path 'STRUM_Log_*.txt'
    
//     shell:
// 	"""
// 	path=$( echo ${path_to_in%/*} )
// 	file=$( echo ${path_to_in##*/} )
//     ./runSTRUM.pl -datadir $path -input_type $file -mutation !{mutation}
//     """
	
// }

// process ROSETTA_DDG_PREMINIMIZATION {
// 	cpus 4
// 	container "${simgDir}/rosetta_23_45_updated_03.sif"

// 	input:
// 	path PDB

// 	output:
// 	path "*.pdb"
// 	path "*.cst"

// 	shell:
// 	"""
// 	minimize_with_cst.mpi.linuxgccrelease -in:file:s ${PDB} -in:file:fullatom -ignore_unrecognized_res -ignore_zero_occupancy false -fa_max_dis 9.0 -database /rosetta/Rosetta/main/database -ddg::harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false > mincst.log
// 	grep 'c-alpha' mincst.log | awk '{print "AtomPair CA "\$8" CA "\$10" HARMONIC "\$12" "\$15""}' > input.cst
// 	"""
// }

// process ROSETTA_DDG {
// 	cpus 4
// 	container "${simgDir}/rosetta_23_45_updated_03.sif"

// 	input:
// 	path min_pdb
// 	path cst
// 	path resf


// 	output:
// 	path "ddg_predictions.out"
// 	path "*_wt.out"
// 	path "*_mt.out"


// 	"""
// 	ddg_monomer.mpi.linuxgccrelease -ignore_zero_occupancy false -in:file:s ${min_pdb} -ddg::mut_file ${resf} -ddg:weight_file soft_rep_design -database /rosetta/Rosetta/main/database -ddg::iterations 50 -ddg::dump_pdbs true -ignore_unrecognized_res -ddg::local_opt_only false -ddg::min_cst true -constraints::cst_file ${cst} -in::file::fullatom -ddg::min true
// 	"""
// }


process MAESTRO {
	cpus 4
	container "${simgDir}/maestro.sif"


	input:
	path effiles_dir
	path council_dir
	path pdb
	val mutation
	val chain
	path xml

	output:
	path "*.csv"

	script:
	mut = mutation[-1]
	mut_index = mutation.indexOf(mut)
	mutation = mutation[0..mut_index-1]
	new_mutation = mutation + '.' + chain + '{' + mut + '}'
	"""
	maestro ${xml} ${pdb} --evalmut=${new_mutation} > maestro_out_${pdb}.csv
	"""

}

process CUTANDMUTATE {
    input:
    path sequence
	val start_pos
	val end_pos
	val mutation
 
    output:
    path 'new_seq.fasta'
 
    """
    mutateAndCut.py -seq "${sequence}" -mutation "${mutation}" -start_position ${start_pos} -end_position ${end_pos} 
    """
}

process CREATEXML {
	input:
	val name
	path sequence
	val pack_round
	val seq_mode
	val skip_unknown_mutant
	val scorefxn
	val start_pos
	val neighbor_dis
	val pack_neighbors
	val weight
	path template

	output:
	path 'thread.xml'

	"""
	sequence_content=\$(cat ${sequence})

	XMLcreate.py -name ${name} -sequence "\${sequence_content}" -start_pos ${start_pos} -pack_neighbors ${pack_neighbors} -neighbor_dis ${neighbor_dis} -scorefxn ${scorefxn} -skip_unknown_mutant ${skip_unknown_mutant} -pack_rounds ${pack_round} -sequence_mode ${seq_mode} -weight ${weight} -template ${template}
	"""
}

process ROSETTA_THREADER {
	cpus 4
	container "${simgDir}/rosetta_23_45_updated_03.sif"

	input:
	path PDB
	path xml_file

	output:
	path "*.pdb"
	path "*.sc"

	script:
	"""
	rosetta_scripts.mpi.linuxgccrelease -s ${PDB} -parser:protocol ${xml_file}
	"""
}

process SPLITPDB {
    input:
    path list_of_structs

    output:
    path "split_chunks/*"

    script:
    """
    # Create a directory to store split chunks
    mkdir -p split_chunks
    
    # Read each line from the input file and split into individual chunks
    cat ${list_of_structs} | awk '{print > "split_chunks/chunk_" NR}'
    """
    
    
}

process ROSETTA_FIXBB {
	cpus 4
	container "${simgDir}/rosetta_23_45_updated_03.sif"
	input:
	path pdb_file 
    path res_file 

	output:
	path "*.pdb"
	path "*.sc"
	path "*.txt"

	script:
	"""
	fixbb.mpi.linuxgccrelease -in:file:s $pdb_file -resfile $res_file -out:suffix _resout > log.txt 
	"""
}


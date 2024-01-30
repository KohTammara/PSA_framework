// r_fbb_PDB = params.list_of_structs
// r_res_fbb = params.res_fbb

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
	val seq_mode
	val pack_round
	val skip_unknown_mutant
	val scorefxn
	val start_pos
	val neighbour_dis
	val pack_neighbours

	output:
	path 'thread.xml'

	script:
	'''
	'''
}

process ROSETTA_THREADER {
	input:
	path path_to_PDB
	val sequence
	val start_pos

	output:
	path 'path_to_PDB_*.pdb'
	path 'score_*.sc'
	/*
	Create script to create xml files for each rosetta threader instance (based on input [read from files might be the answer for both amny proteins for single template and mmany templates])
	edit commands accordingly, may need conditionals and checks for list of templates for a single seq, list of seq for single template, and combinations.
	*/
	script:
	"""
	rosetta_scripts.mpi.linuxgccrelease -s 1ubq.pdb -parser:protocol nothing.xml
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
	container "${simgDir}/rosetta_23_45_reduced.sif"
	
	input:
	path pdb_file 
    path res_file 

	output:
	path "*.pdb"
	path "*.sc"
	path "*.txt"
	/*
	Create script to create xml files for each rosetta threader instance (based on input [read from files might be the answer])
	*/
	script:
	"""
	fixbb.mpi.linuxgccrelease -in:file:s $pdb_file -resfile $res_file -out:suffix _resout > log.txt 
	"""
}


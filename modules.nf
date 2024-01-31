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
	container "${simgDir}/rosetta_23_45_reduced.sif"

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


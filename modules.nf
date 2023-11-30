process STRUM {
	/*
	Adapt to read from file for multiple instances
	*/
	input:
	val mutation
	val input_type
	path path_to_in
	
	output:
    path 'STRUM_Log_*.txt'
    
    shell:
	"""
	path=$( echo ${path_to_in%/*} )
	file=$( echo ${path_to_in##*/} )
    ./runSTRUM.pl -datadir $path -input_type $file -mutation !{mutation}
    """
	
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

process ROSETTA_FIXBB {
	input:
	path path_to_PDB
	path path_to_res

	output:
	path 'path_to_PDB_*.pdb'
	path 'score_*.sc'
	path 'log_*.txt'
	/*
	Create script to create xml files for each rosetta threader instance (based on input [read from files might be the answer])
	*/
	script:
	"""
	fixbb.mpi.linuxgccrelease -in:file:s $path_to_PDB -resfile $path_to_res -out:suffix _resout > log_*.txt 
	"""
}

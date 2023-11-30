process STRUM {
	input:
	val mutation
	val input_type
	path path_to_in
	path path_to_out
	
	output:
    path 'STRUM_Log_*'
    
    script:
    """
    ./runSTRUM.pl -datadir path_to_in -input_type [some variable holding input] -mutation mutation
    """
	
}

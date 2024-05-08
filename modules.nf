//create a collective gromacs process or see if it can be created otherwise split it at points where collective breaks
//process for native and a process for mutant
process GROMACS_MT_FBB {
	publishDir '/GROMACS_FBB_mutant', mode: 'copy', overwrite: false
	container "${simgDir}/gromacs2023_2_mpi_charmm36m.sif"

	input:
	path pdb_file
	path ions_mdp
	path em_mdp
	path nvt_mdp
	path npt_mdp
	path md_mdp
	val G1
	val G2
	val G3
	val G4
	val G5
	val G6
	val G7
	val G8
	val G9

	output:
	path '*_processed.gro'
	path '*.top'
	path '*_box.gro'
	path '*_solve.gro'
	path '*_ions.tpr'
	path '*_em.tpr'
	path '*_nvt.tpr'

	script:
	"""
	grep -v HOH ${pdb_file} > ${pdb_file}_clean.pdb
	echo ${G1} | gmx_mpi pdb2gmx -f ${pdb_file}_clean.pdb -o ${pdb_file}_processed.gro -water spce -ignh
	gmx_mpi editconf -f ${pdb_file}_processed.gro -o ${pdb_file}_box.gro -c -d 1.0 -bt dodecahedron
	gmx_mpi solvate -cp ${pdb_file}_box.gro -cs spc216.gro -o  ${pdb_file}_solve.gro -p topol.top
	gmx_mpi grompp -f ${ions_mdp} -c ${pdb_file}_solve.gro -p topol.top -o ${pdb_file}_ions.tpr -maxwarn 1
	echo ${G2} |gmx_mpi genion -s ${pdb_file}_ions.tpr -o ${pdb_file}_solve_ions.gro -p topol.top -pname SOD -nname CLA -neutral
	gmx_mpi grompp -f ${em_mdp} -c ${pdb_file}_solve_ions.gro -p topol.top -o ${pdb_file}_em.tpr
	gmx_mpi mdrun -v -deffnm ${pdb_file}_em
	echo ${G3}|gmx_mpi energy -f ${pdb_file}_em.edr -o ${pdb_file}_potential.xvg
	gmx_mpi grompp -f ${nvt_mdp} -c ${pdb_file}_em.gro -r ${pdb_file}_em.gro -p topol.top -o ${pdb_file}_nvt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_nvt
	echo ${G4}|gmx_mpi energy -f ${pdb_file}_nvt.edr -o ${pdb_file}_temperature.xvg
	gmx_mpi grompp -f ${npt_mdp} -c ${pdb_file}_nvt.gro -r ${pdb_file}_nvt.gro -t ${pdb_file}_nvt.cpt -p topol.top -o ${pdb_file}_npt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_npt
	echo ${G5}|gmx_mpi energy -f ${pdb_file}_npt.edr -o ${pdb_file}_density.xvg
	gmx_mpi grompp -f ${md_mdp} -c ${pdb_file}_npt.gro -t ${pdb_file}_npt.cpt -p topol.top -o ${pdb_file}_md_0_1.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_md_0_1
	echo ${G6} |gmx_mpi trjconv -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1.xtc -o ${pdb_file}_md_0_1_noPBC.xtc -pbc mol -center
	echo ${G7} |gmx_mpi rms -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd.xvg -tu ns
	echo ${G8} |gmx_mpi rms -s ${pdb_file}_em.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd_xtal.xvg -tu ns
	echo ${G9} |gmx_mpi gyrate -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_gyrate.xvg
	"""


}

process GROMACS_MT_THREADER {
	publishDir '${params.pub_dir}/GROMACS_THREADER_mutant', mode: 'copy', overwrite: false
	container "${simgDir}/gromacs2023_2_mpi_charmm36m.sif"

	input:
	path pdb_file
	path ions_mdp
	path em_mdp
	path nvt_mdp
	path npt_mdp
	path md_mdp
	val G1
	val G2
	val G3
	val G4
	val G5
	val G6
	val G7
	val G8
	val G9

	output:
	path '*_processed.gro'
	path '*.top'
	path '*_box.gro'
	path '*_solve.gro'
	path '*_ions.tpr'
	path '*_em.tpr'
	path '*_nvt.tpr'

	script:
	"""
	grep -v HOH ${pdb_file} > ${pdb_file}_clean.pdb
	echo ${G1} | gmx_mpi pdb2gmx -f ${pdb_file}_clean.pdb -o ${pdb_file}_processed.gro -water spce -ignh
	gmx_mpi editconf -f ${pdb_file}_processed.gro -o ${pdb_file}_box.gro -c -d 1.0 -bt dodecahedron
	gmx_mpi solvate -cp ${pdb_file}_box.gro -cs spc216.gro -o  ${pdb_file}_solve.gro -p topol.top
	gmx_mpi grompp -f ${ions_mdp} -c ${pdb_file}_solve.gro -p topol.top -o ${pdb_file}_ions.tpr -maxwarn 1
	echo ${G2} |gmx_mpi genion -s ${pdb_file}_ions.tpr -o ${pdb_file}_solve_ions.gro -p topol.top -pname SOD -nname CLA -neutral
	gmx_mpi grompp -f ${em_mdp} -c ${pdb_file}_solve_ions.gro -p topol.top -o ${pdb_file}_em.tpr
	gmx_mpi mdrun -v -deffnm ${pdb_file}_em
	echo ${G3}|gmx_mpi energy -f ${pdb_file}_em.edr -o ${pdb_file}_potential.xvg
	gmx_mpi grompp -f ${nvt_mdp} -c ${pdb_file}_em.gro -r ${pdb_file}_em.gro -p topol.top -o ${pdb_file}_nvt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_nvt
	echo ${G4}|gmx_mpi energy -f ${pdb_file}_nvt.edr -o ${pdb_file}_temperature.xvg
	gmx_mpi grompp -f ${npt_mdp} -c ${pdb_file}_nvt.gro -r ${pdb_file}_nvt.gro -t ${pdb_file}_nvt.cpt -p topol.top -o ${pdb_file}_npt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_npt
	echo ${G5}|gmx_mpi energy -f ${pdb_file}_npt.edr -o ${pdb_file}_density.xvg
	gmx_mpi grompp -f ${md_mdp} -c ${pdb_file}_npt.gro -t ${pdb_file}_npt.cpt -p topol.top -o ${pdb_file}_md_0_1.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_md_0_1
	echo ${G6} |gmx_mpi trjconv -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1.xtc -o ${pdb_file}_md_0_1_noPBC.xtc -pbc mol -center
	echo ${G7} |gmx_mpi rms -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd.xvg -tu ns
	echo ${G8} |gmx_mpi rms -s ${pdb_file}_em.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd_xtal.xvg -tu ns
	echo ${G9} |gmx_mpi gyrate -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_gyrate.xvg
	"""


}

process GROMACS_WT {
	publishDir "${params.pub_dir}/GROMACS_wildtype", mode: 'copy', overwrite: false
	container "${simgDir}/gromacs2023_2_mpi_charmm36m.sif"

	input:
	path pdb_file
	path ions_mdp
	path em_mdp
	path nvt_mdp
	path npt_mdp
	path md_mdp
	val G1
	val G2
	val G3
	val G4
	val G5
	val G6
	val G7
	val G8
	val G9

	output:
	path '*_processed.gro'
	path '*.top'
	path '*_box.gro'
	path '*_solve.gro'
	path '*_ions.tpr'
	path '*_em.tpr'
	path '*_nvt.tpr'

	script:
	"""
	grep -v HOH ${pdb_file} > ${pdb_file}_clean.pdb
	echo ${G1} | gmx_mpi pdb2gmx -f ${pdb_file}_clean.pdb -o ${pdb_file}_processed.gro -water spce -ignh
	gmx_mpi editconf -f ${pdb_file}_processed.gro -o ${pdb_file}_box.gro -c -d 1.0 -bt dodecahedron
	gmx_mpi solvate -cp ${pdb_file}_box.gro -cs spc216.gro -o  ${pdb_file}_solve.gro -p topol.top
	gmx_mpi grompp -f ${ions_mdp} -c ${pdb_file}_solve.gro -p topol.top -o ${pdb_file}_ions.tpr -maxwarn 1
	echo ${G2} |gmx_mpi genion -s ${pdb_file}_ions.tpr -o ${pdb_file}_solve_ions.gro -p topol.top -pname SOD -nname CLA -neutral
	gmx_mpi grompp -f ${em_mdp} -c ${pdb_file}_solve_ions.gro -p topol.top -o ${pdb_file}_em.tpr
	gmx_mpi mdrun -v -deffnm ${pdb_file}_em
	echo ${G3}|gmx_mpi energy -f ${pdb_file}_em.edr -o ${pdb_file}_potential.xvg
	gmx_mpi grompp -f ${nvt_mdp} -c ${pdb_file}_em.gro -r ${pdb_file}_em.gro -p topol.top -o ${pdb_file}_nvt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_nvt
	echo ${G4}|gmx_mpi energy -f ${pdb_file}_nvt.edr -o ${pdb_file}_temperature.xvg
	gmx_mpi grompp -f ${npt_mdp} -c ${pdb_file}_nvt.gro -r ${pdb_file}_nvt.gro -t ${pdb_file}_nvt.cpt -p topol.top -o ${pdb_file}_npt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_npt
	echo ${G5}|gmx_mpi energy -f ${pdb_file}_npt.edr -o ${pdb_file}_density.xvg
	gmx_mpi grompp -f ${md_mdp} -c ${pdb_file}_npt.gro -t ${pdb_file}_npt.cpt -p topol.top -o ${pdb_file}_md_0_1.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_md_0_1
	echo ${G6} |gmx_mpi trjconv -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1.xtc -o ${pdb_file}_md_0_1_noPBC.xtc -pbc mol -center
	echo ${G7} |gmx_mpi rms -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd.xvg -tu ns
	echo ${G8} |gmx_mpi rms -s ${pdb_file}_em.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd_xtal.xvg -tu ns
	echo ${G9} |gmx_mpi gyrate -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_gyrate.xvg
	"""


}


process MAESTRO_XML {
	input:
	path effiles_dir
	path council_dir
	val path_to_pdb
	val prefix
	val postfix
	val to_lower
	val bu

	output:
	path "config.xml"


	"""
	xml_maestro.py -pdb_path "${path_to_pdb}" -prefix "${prefix}" -postfix "${postfix}" -tolower "${to_lower}" -bu "${bu}" > config.xml
	"""
}

process MAESTRO {
	publishDir "${params.pub_dir}/Maestro/out", mode: 'copy', overwrite: true, enabled: true
	cpus 6
	container "${simgDir}/maestro.sif"


	input:
	path effiles_dir
	path council_dir
	path mutation
	path xml

	output:
	path "*.csv"

	script:
	"""
	maestro ${xml} --evalmutlist=${mutation} --energy > maestro_out.csv
	"""

}

process CUTANDMUTATE {
    input:
    path sequence
	path mutation_info
 
    output:
    path 'new_seq.fasta'
	val 'start_position.txt'
 
    """
    mutateAndCut.py -seq "${sequence}" -mutation "${mutation_info}" > start_position.txt
    """
}

process CREATEXML {
	conda 'env.yaml'

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
	publishDir "${params.pub_dir}/Rosetta_threader/${mutation}", mode: 'copy', overwrite: false
	container "${simgDir}/rosetta_23_45_updated_03.sif"

	input:
	path PDB
	path xml_file
	val mutation

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
	publishDir "${params.pub_dir}/Rosetta_FixBB/${res_file}", mode: 'copy', overwrite: false
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

process PMX_PREP {
	publishDir "${params.pub_dir}/pmx_preparation/${params.pdb}/${params.res_number}_${params.mutant_res}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx.sif"

	input:
	path pdb
	val res_number
	val res_mutant
	val forcefield

	output:
	path "mutant.pdb"
	path "newtop.top"

	script:
	//redo as CLI and not py script, api weird.
	"""
	python3 ${params.bin_dir}/pmx_prep.py -pdb '${pdb}' -res_num ${res_number} -res_mut '${res_mutant}' -ff '${forcefield}'
	"""
}


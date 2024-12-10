/*
All processes used are defined here.
Processes are called in main.nf for use independantly or within a workflow.
Directories where results from processes are published can be changed by changing the path of publishDir in params.json.
In the event that differing containers are used, the name of the sif file should be changed.
	Containers are stored in the image directory of this project.
	The name will following `container` directive. 
Descriptions of each process given after process definition. 
*/

process GROMACS_MT_PMX {
/*
	Classical MD simulation of mutant folded structure (Structure mutated with pmx). 
	Input : 
		MDP files for energy minisation (em_mdp), ion generation (ions_mdp),
		NVT (nvt_em), NPT (npt_em), production MD (md_mdp,
		the structure file (pdb_file), and the various selections from input prompts 
		(details of the selection options can be found in GROMACS_processes_info.txt in the dirrectory containing this file)
	Output :
		All processed output files such as .gro/.tpr/.top files created from gmx_mpi commands.
		All xvg files (potential/temperature/density/rmsd/gyration)
*/
	publishDir "${params.pub_dir}/GROMACS_pmx_mutant", mode: 'copy', overwrite: false
	container "${simgDir}/gromacs2023_2_mpi_charmm36m.sif"
	tag "${pdb_file}"
	memory '8GB'
	cpus '8'

	input:
	path pdb_file
	path topol
	path itp
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
	val G10

	output:
	path '*.gro'
	path '*.tpr'
	path "${topol}"
	path '*_box.gro'
	path '*_solve.gro'
	path '*_ions.tpr'
	path '*_em.tpr'
	path '*_nvt.tpr'
	path '*.xvg'

	script:
	"""
	gmx_mpi grompp -f ${em_mdp} -c ${pdb_file} -p ${topol} -o ${pdb_file}_em.tpr
	gmx_mpi mdrun -v -deffnm ${pdb_file}_em
	echo ${G3}|gmx_mpi energy -f ${pdb_file}_em.edr -o ${pdb_file}_potential.xvg
	gmx_mpi grompp -f ${nvt_mdp} -c ${pdb_file}_em.gro -r ${pdb_file}_em.gro -p ${topol} -o ${pdb_file}_nvt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_nvt
	echo ${G4}|gmx_mpi energy -f ${pdb_file}_nvt.edr -o ${pdb_file}_temperature.xvg
	gmx_mpi grompp -f ${npt_mdp} -c ${pdb_file}_nvt.gro -r ${pdb_file}_nvt.gro -t ${pdb_file}_nvt.cpt -p ${topol} -o ${pdb_file}_npt.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_npt
	echo ${G5} | gmx_mpi energy -f ${pdb_file}_npt.edr -o ${pdb_file}_pressure.xvg
	echo ${G6}|gmx_mpi energy -f ${pdb_file}_npt.edr -o ${pdb_file}_density.xvg
	gmx_mpi grompp -f ${md_mdp} -c ${pdb_file}_npt.gro -t ${pdb_file}_npt.cpt -p ${topol} -o ${pdb_file}_md_0_1.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_md_0_1
	echo ${G7} |gmx_mpi trjconv -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1.xtc -o ${pdb_file}_md_0_1_noPBC.xtc -pbc mol -center
	echo ${G8} |gmx_mpi rms -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd.xvg -tu ns
	echo ${G9} |gmx_mpi rms -s ${pdb_file}_em.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd_xtal.xvg -tu ns
	echo ${G10} |gmx_mpi gyrate -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_gyrate.xvg
	"""


}

process STANDARD_MD {
/*
	Classical MD simulation of mutant folded structure (Structure mutated with Rosetta FBB). 
	Input : 
		MDP files for energy minisation (em_mdp), ion generation (ions_mdp),
		NVT (nvt_em), NPT (npt_em), production MD (md_mdp,
		the structure file (pdb_file), and the various selections from input prompts 
		(details of the selection options can be found in GROMACS_processes_info.txt in the dirrectory containing this file)
	Output :
		All processed output files such as .gro/.tpr/.top files created from gmx_mpi commands.
		All xvg files (potential/temperature/density/rmsd/gyration)
*/
	publishDir "${params.pub_dir}/GROMACS_STD_mutant/", mode: 'copy', overwrite: false
	container "${simgDir}/gromacs2023_2_mpi_charmm36m.sif"
	tag "${pdb_file}"

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
	val G10

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
	echo ${G5} | gmx_mpi energy -f ${pdb_file}_npt.edr -o ${pdb_file}_pressure.xvg
	echo ${G6}|gmx_mpi energy -f ${pdb_file}_npt.edr -o ${pdb_file}_density.xvg
	gmx_mpi grompp -f ${md_mdp} -c ${pdb_file}_npt.gro -t ${pdb_file}_npt.cpt -p topol.top -o ${pdb_file}_md_0_1.tpr
	gmx_mpi mdrun -deffnm ${pdb_file}_md_0_1
	echo ${G7} |gmx_mpi trjconv -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1.xtc -o ${pdb_file}_md_0_1_noPBC.xtc -pbc mol -center
	echo ${G8} |gmx_mpi rms -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd.xvg -tu ns
	echo ${G9} |gmx_mpi rms -s ${pdb_file}_em.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_rmsd_xtal.xvg -tu ns
	echo ${G10} |gmx_mpi gyrate -s ${pdb_file}_md_0_1.tpr -f ${pdb_file}_md_0_1_noPBC.xtc -o ${pdb_file}_gyrate.xvg
	"""


}


process MAESTRO_XML {
/*
	Process to create the Maestro XML configuration file through the use of a python script.
	Input : 
		Path to effiles and council directories required for maestro execution.
		Path to the PDB file, prefix for output, postfix for output, switch for lowercase,
		bu is a switch for if the directory contains biological assemblies.
	Output :
		Maestro configuration file.
*/
	publishDir "${params.pub_dir}/Maestro_xml", mode: 'copy', overwrite: false
	memory '100MB'
	cpus '4'

	input:
	path effiles_dir
	path council_dir
	val path_to_pdb
	val prefix
	val postfix
	val to_lower
	val to_upper
	val bu

	output:
	path "config.xml"


	"""
	xml_maestro.py -pdb_path "${path_to_pdb}" -prefix "${prefix}" -postfix "${postfix}" -tolower "${to_lower}" -toupper "${to_upper}" -bu "${bu}" > config.xml
	"""
}

process MAESTRO {
/*
	Evaluate mutants with Maestro which provides free energy values in the form of a CSV file as output.
	Input : 
		XML configuration file that is obtained from process MAESTRO_XML.
		Path to effiles and council directories required for Maestro execution.
		Path to file containing mutations to be evaluated using Maestro.
	Output :
		CSV containing mutant free energy and scores in format:
		structure<TAB>seqlength<TAB>mutation<TAB>score<TAB>delta_score<TAB>ddG<TAB>ddG_confidence
*/
	publishDir "${params.pub_dir}/Maestro/out", mode: 'copy', overwrite: false
	cpus '4'
	container "${simgDir}/maestro.sif"
	tag "${mutation}"


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
/*
	Cut and insert mutation according to specified range and mutant given through the use of a Python script
	Input : 
		Path to file containing the sequence(Fasta file).
		Path to file containing the mutation information.
	Output :
		Fasta file containing new cut and mutated sequence.
		Start position of sequence taken from mutation information.
*/
	publishDir "${params.pub_dir}/cut_and_mutate_sequence", mode: 'copy', overwrite: false
	memory '100MB'
	cpus '4'
	tag "${mutation_info}"

    input:
    path sequence
	path mutation_info
 
    output:
    path 'new_seq.fasta'
	val '1'
 
    """
    mutateAndCut.py -seq "${sequence}" -mutation "${mutation_info}" 
    """
}

process CREATEXML {
/*
	Makes use of a python script to create the XML file required as input for the Rosetta SImpleThreadingMover
	Input : 
		Name of threading process.
		Path to sequence to be used in the threading process (Fasta File).
		Pack_round, seq_mode, skip_unknown_mutant, scorefxn, start_pos, neighbor_dis, pack_neighbors, and weight are variables 
		in the created XML file and description of these can be found in rosetta commons under SimpleThradingMover.
		Path to template structure that the sequence is threaded onto.
	Output :
		XML file to execute the Rosetta SimpleThreadingMover 
*/
	publishDir "${params.pub_dir}/Rosetta_threader_xml", mode: 'copy', overwrite: false
	memory '100MB'
	cpus '4'
	conda 'env.yaml'
	// beforeScript "/apps2/mambaforge/bin/conda env create --file /home/21616019/git_projects/v5/PSA_framework/env.yaml"
	// beforeScript "/apps2/mambaforge/bin/conda init"
	// beforeScript "/apps2/mambaforge/bin/conda activate my-env-1"
	// beforeScript 'alias python=/apps2/mambaforge/envs/my-env-1/bin/python3'
	// afterScript "/apps2/mambaforge/bin/conda deactivate"


	input:
	val name
	path sequence
	val pack_round
	val seq_mode
	val skip_unknown_mutant
	val scorefxn
	path start_pos
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
/*
	Process to execute the Rosetta SimpleThreadingMover
	Input : 
		Path to PDB file to be used as a template structure
		XML file created from process CREATEXML
		Path to file containing mutation information (Name of file is used as a tag for the process)
	Output :
		Path to PDB file containing the threaded structure.
		Path to score file for threaded structure
		Name of mutation taken from mutation file
*/
	publishDir "${params.pub_dir}/Rosetta_threader/${mutation.baseName}", mode: 'copy', overwrite: false
	container "${simgDir}/rosetta_23_45_updated_03.sif"
	memory '1GB'
	cpus '8'
	tag "${mutation.baseName}"

	input:
	path PDB
	path xml_file
	path mutation

	output:
	path "*.pdb"
	path "*.sc"
	val "${mutation.baseName}"

	script:
	"""
	rosetta_scripts.mpi.linuxgccrelease -s ${PDB} -parser:protocol ${xml_file}
	"""
}


process ROSETTA_FIXBB {
/*
	This process executes the rosetta FixedBackBone application
	Input : 
		Path to PDB file of structure to be mutated.
		Residue file used by the rosetta FBB application (Describes mutations and state of backbone)
	Output :
		PDB file of mutated structure from Rosetta FBB
		Score file from Rosetta FBB application
		Text file containing a log of the Rosetta FBB process
		Name of the residue file which is used to take the process and tag following processes
*/
	publishDir "${params.pub_dir}/Rosetta_FixBB/${res_file.baseName}", mode: 'copy', overwrite: false
	container "${simgDir}/rosetta_23_45_updated_03.sif"
	tag "${res_file.baseName}"
	memory '1GB'
	cpus '8'
	
	input:
	path pdb_file
	path res_file


	output:
	path "*.pdb" 
	path "*.sc" 
	path "*.txt"
	val "${res_file.baseName}"

	script:
	"""
	fixbb.mpi.linuxgccrelease -in:file:s $pdb_file -resfile $res_file -out:suffix _resout > log.txt 
	"""
}

process PMX_PREP_MUTANT {
/*
	This process selects the appropriate forcefield dependant on user input and makes use of PMX to mutate structures for ion generation,
	placement of the protein in a box and solvation.
	Input:
		Path to pmx mutant PBD file
		Tuple containing the index of residue to be mutated and mutant residue
		forcefield to be used
		Path to mdp file for ion generation (genion)
	Output:
		Mutant pdb file
		updated topology file
		topology file generated by pdb2gmx
		pdb file of structure in a box
		Solvated pdb 
		tpr file from genion
		ions.pdb is pdb file after neutralisation
		posre.itp position restraints file
		value containing index and mutant for naming further in the workflow

*/
	publishDir "${params.pub_dir}/pmx/preparation/${index}_${mutant}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	tag "${index}_${mutant}"
	memory '1GB'
	cpus '8'

	input:
	path pdb
	tuple val(index), val(mutant)
	val forcefield
	path genion

	output:
	path "mutant.pdb"
	path "newtop.top"
	path "conf.pdb"
	path "box.pdb"
	path "water.pdb"
	path "genion.tpr"
	path "${index}_${mutant}_ions.pdb"
	path "posre.itp"
	val "${index}_${mutant}"

	script:
	def forcefield_number
	switch(forcefield) {
		case 'amber14sbmut':
			forcefield_number = '1'
			break
		case 'amber99sb-star-ildn-bsc1-mut':
			forcefield_number = '2'
			break
		case 'amber99sb-star-ildn-dna-mut':
			forcefield_number = '3'
			break
		case 'amber99sb-star-ildn-mut':
			forcefield_number = '4'
			break
		case 'charmm22star-mut ':
			forcefield_number = '5'
			break
		case 'charmm36m-mut':
			forcefield_number = '6'
			break
	}
	index = index-124
	"""
	echo -e "${forcefield_number}\n${index}\n${mutant}\nn\n" | pmx mutate -f ${pdb} -o mutant.pdb
	gmx_mpi pdb2gmx -f mutant.pdb -o conf.pdb -p topol.top -ff ${forcefield} -water tip3p
	sed -i 's#/usr/local/lib/python3.10/dist-packages/pmx/data/mutff/##g' topol.top
	pmx gentop -p topol.top -o newtop.top
	gmx_mpi editconf -f conf.pdb -o box.pdb -bt dodecahedron -d 1.0
	gmx_mpi solvate -cp box -cs spc216 -p newtop -o water.pdb
	gmx_mpi grompp -f ${genion} -c water.pdb -p newtop.top -o genion.tpr
	echo "SOL" | gmx_mpi genion -s genion.tpr -p newtop.top -neutral -conc 0.15 -o ${index}_${mutant}_ions.pdb
	"""
}

process GRO_PREP_MUTANT {
/*
	This process makes use of Rosetta mutated structures and prepares them for a classical simulation with ion generation,
	placement of the protein in a box and solvation.
	Input:
		Path to pmx mutant PBD file
		Tuple containing the index of residue to be mutated and mutant residue
		forcefield to be used
		Path to mdp file for ion generation (genion)
	Output:
		updated topology file
		topology file generated by pdb2gmx
		pdb file of structure in a box
		Solvated pdb 
		tpr file from genion
		ions.pdb is pdb file after neutralisation
		posre.itp position restraints file
		value containing index and mutant for naming further in the workflow

*/
	publishDir "${params.pub_dir}/gro_preparation_rosetta/${mutation}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	tag "${mutation}"


	input:
	path pdb
	val mutation
	val ff_name
	path genion

	output:
	path "topol.top" 
	path "conf.pdb"
	path "box.pdb"
	path "water.pdb"
	path "genion.tpr"
	path "ions.pdb"
	path "posre.itp"
	val "${mutation}"

	script:
	"""
	gmx_mpi pdb2gmx -f ${pdb} -o conf.pdb -p topol.top -ff ${ff_name} -water tip3p -ignh
	gmx_mpi editconf -f conf.pdb -o box.pdb -bt dodecahedron -d 1.0
	gmx_mpi solvate -cp box -cs spc216 -p topol -o water.pdb
	gmx_mpi grompp -f ${genion} -c water.pdb -p topol.top -o genion.tpr
	echo "SOL" | gmx_mpi genion -s genion.tpr -p topol.top -neutral -conc 0.15 -o ions.pdb
	"""

}

process GRO_PREP_TRIPEPTIDE {
/*
	This process PMX tripeptide structures for ion generation,
	placement of the protein in a box and solvation.
	Input:
		Path to tuple containing the tripeptide topology file, position restraints file (posre.itp), coordinate file (itp), and name of tripeptide.
		Path to mdp file for ion generation (genion)
		forcefield to be used (not exactly necessary in this case as ff is specified in topology)
	Output:
		updated topology file
		pdb file for tripeptide
		pdb file of structure in a box
		Solvated pdb 
		tpr file from genion
		ions.pdb is pdb file after neutralisation
		posre.itp position restraints file
		coordinate file(itp)
		name of tripeptide

*/
	publishDir "${params.pub_dir}/gro_preparation/tripeptide/${name}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	tag "${name}"
	memory '100MB'
	cpus '8'
	time '5m'

	input:
	tuple path(posre_itp),	path(top), path(pdb), val(name), path(itp)
	path genion
	val ff_name
	

	output:
	path "${top}"
	path "${pdb}"
	path "box.pdb"
	path "water.pdb"
	path "genion.tpr"
	path "ions.pdb"
	path "${posre_itp}"
	path "${itp}"
	val name

	script:
	"""
	sed -i 's/charmm36mut.ff/charmm36m-mut.ff/g' ${top}
    sed -i 's/tips3p.itp/tip3p.itp/g' ${top}
	gmx_mpi editconf -f ${pdb} -o box.pdb -bt dodecahedron -d 1.0
	gmx_mpi solvate -cp box -cs spc216 -p ${top} -o water.pdb
	gmx_mpi grompp -f ${genion} -c water.pdb -p ${top} -o genion.tpr -maxwarn 1
	echo "SOL" | gmx_mpi genion -s genion.tpr -o ions.pdb -p ${top} -pname NA -nname CL -neutral
	"""
}

process READ_TRIPEPTIDE_FILES {
/*
	This process obtains the files of a tripeptide gathered from the pmx tripeptide database. These files incude  the topology file, itp file and position restraints (itp) file.
	Input:
		Path to directory containing the folders for tripeptides.
	Output:
		Tuple containing the 3 files for each tripeptide. 
*/
	tag "${tripep_dir.baseName}"
	memory '100MB'
	cpus '8'
	time '5m'

	input:
	path tripep_dir

	output:
	tuple path("${tripep_dir}/posre*.itp"),	path("${tripep_dir}/*.top"), path("${tripep_dir}/*.pdb"), val("$tripep_dir.baseName"), path("${tripep_dir}/${tripep_dir.baseName}.itp")

	script:
	"""
    itp_file1=''
    itp_file2=''
    top_file=''
    pdb_file=''
    
    for file in ${tripep_dir}/*; do
        case \$file in
            ${tripep_dir}/posre*.itp) itp_file1=\$file ;;
            ${tripep_dir}/*.itp) [[ -z "\$itp_file2" && "\$file" != "\$itp_file1" ]] && itp_file2=\$file ;;
            ${tripep_dir}/*.top) top_file=\$file ;;
            ${tripep_dir}/*.pdb) pdb_file=\$file ;;
        esac
    done

	# Create symbolic links for the output files
    ln -s "\$itp_file1" itp_file1
    ln -s "\$itp_file2" itp_file2
    ln -s "\$top_file" top_file
    ln -s "\$pdb_file" pdb_file
	""" 



}

process GRO_EQUILIBRIUM {
/*
	Process executing energy minimisation, npt and equilibration of protein for free energy simulation. This process is for the folded structures (mutants).
	Input:
		Ions_pdb generated from preperation (genion neutralisation)
		topol topology file
		posre_itp position restraits file
		energy minimisation mdp file
		equilibration mdp file
		npt mdp file
		name of mutant
	Output:
		equilibration trajectory file
		equilibration run input file
		energy minimisation run input file
		npt run input file
		name of mutant
*/
	publishDir "${params.pub_dir}/gro_preparation/equilibrium/${name}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	tag "${name}"
	memory '12GB'
	cpus '12'
	time '3d'

	input:
	path ions_pdb
	path topol
	path posre_itp
	path enmin
	path equil
	path npt
	val name

	output:
	path "equil.trr"
	path "equil.tpr"
	path "enmin.tpr"
	path "npt.tpr"
	val "${name}"
	path "${topol}"

	script:
	"""
	gmx_mpi grompp -f ${enmin} -c ${ions_pdb} -p ${topol} -o enmin.tpr -maxwarn 1
	gmx_mpi mdrun -s enmin.tpr -deffnm enmin -v
	gmx_mpi grompp -f ${npt} -c enmin.gro -r enmin.gro -p ${topol} -o npt.tpr -maxwarn 2
	gmx_mpi mdrun -s npt.tpr -deffnm npt -v
	gmx_mpi grompp -f ${equil} -c npt.gro -p ${topol} -o equil.tpr -maxwarn 1
	gmx_mpi mdrun -s equil.tpr -deffnm equil -v
	"""
}

process GRO_EQUILIBRIUM_UNFOLDED {
/*
	Process executing energy minimisation, npt and equilibration of protein for free energy simulation. This process is for the unfolded structures (tripeptide).
	Input:
		Ions_pdb generated from preperation (genion neutralisation)
		topol topology file
		posre_itp position restraits file
		energy minimisation mdp file
		equilibration mdp file
		npt mdp file
		name of mutant
	Output:
		equilibration trajectory file
		equilibration run input file
		energy minimisation run input file
		npt run input file
		name of mutant
		topol topology file
		position restraint ipt file
		itp file
*/
	publishDir "${params.pub_dir}/gro_preparation/equilibrium_unfolded/${name}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	tag "${name}"
	memory '32GB'
	cpus '12'
	time '3d'

	input:
	path ions_pdb
	path topol
	path posre_itp
	path itp
	path enmin
	path equil
	path npt
	val name

	output:
	path "equil.trr"
	path "equil.tpr"
	path "enmin.tpr"
	path "npt.tpr"
	val "${name}"
	path "${topol}"
	path "${posre_itp}"
	path "${itp}"

	script:
	"""
	gmx_mpi grompp -f ${enmin} -c ${ions_pdb} -p ${topol} -o enmin.tpr -maxwarn 2
	gmx_mpi mdrun -s enmin.tpr -deffnm enmin -v
	gmx_mpi grompp -f ${npt} -c enmin.gro -r enmin.gro -p ${topol} -o npt.tpr -maxwarn 3
	gmx_mpi mdrun -s npt.tpr -deffnm npt -v
	gmx_mpi grompp -f ${equil} -c npt.gro -p ${topol} -o equil.tpr -maxwarn 3
	gmx_mpi mdrun -s equil.tpr -deffnm equil -v
	"""
}

process GRO_NON_EQUILIBRIUM_UNFOLDED {
/*
	Process executing non-equilibrium simulation of protein for free energy simulations. This process is for the unfolded structures (tripeptides).
	Input:
		equilibration trajectory file
		equilibration run input file
		topol topology file
		posre_itp position restraits file
		itp file
		non_equil mdp file
		name of mutant
		type of structure
	Output:
		dgdl xvg files
		name of structure
*/
	publishDir "${params.pub_dir}/gro_preparation/non_equilibrium/${name}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	scratch 'scratch'
	stageInMode 'copy'
	tag "${name}"
	memory '12GB'
	cpus '12'
	time '3d'

	input:
	path equil_trr
	path equil_tpr
	path topol
	path posre_itp
	path itp
	path non_equil
	val name
	val type

	output:
	path "dgdl*", emit: dgdlFiles
	val "${name}"

	script:
	//extract 50 snapshots from the 5ns equilibrium sim (1 per 100ps starting at 100ps)
	"""
	echo "System" | gmx_mpi trjconv -f ${equil_trr} -s ${equil_tpr} -sep -b 100 -o frame_.gro
	for i in \$( seq 0 49 ); do
		n=\$((i+1));
		mkdir frame\$n;
		mv frame_\$i.gro frame\$n/frame.gro;
	done

	for i in \$( seq 1 50 ); do
		cd frame\$i;
		gmx_mpi grompp -f ../${non_equil} -c frame.gro -p ../${topol} -o nonequil.tpr -maxwarn 1;
		gmx_mpi mdrun -s nonequil.tpr -deffnm nonequil -dhdl dgdl_${type}\$i.xvg -v;
		mv dgdl_${type}\$i.xvg ../;
		cd ../;
	done
	"""

}

process GRO_NON_EQUILIBRIUM {
/*
	Process executing non-equilibrium simulation of protein for free energy simulations. This process is for the folded structures (mutant).
	Input:
		equilibration trajectory file
		equilibration run input file
		topol topology file
		non_equil mdp file
		name of mutant
	Output:
		dgdl xvg files
		name of structure
*/
	publishDir "${params.pub_dir}/gro_preparation/non_equilibrium/${name}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	scratch 'scratch'
	stageInMode 'copy'
	tag "${name}"
	memory '12GB'
	cpus '12'
	time '3d'

	input:
	path equil_trr
	path equil_tpr
	path topol
	path non_equil
	val name
	val type

	output:
	path "dgdl*", emit: dgdlFiles
	val "${name}"

	script:
	//extract 50 snapshots from the 5ns equilibrium sim (1 per 100ps starting at 100ps)
	"""
	echo "System" | gmx_mpi trjconv -f ${equil_trr} -s ${equil_tpr} -sep -b 100 -o frame_.gro
	for i in \$( seq 0 49 ); do
		n=\$((i+1));
		mkdir frame\$n;
		mv frame_\$i.gro frame\$n/frame.gro;
	done

	for i in \$( seq 1 50 ); do
		cd frame\$i;
		gmx_mpi grompp -f ../${non_equil} -c frame.gro -p ../${topol} -o nonequil.tpr -maxwarn 1;
		gmx_mpi mdrun -s nonequil.tpr -deffnm nonequil -dhdl dgdl_${type}\$i.xvg -v;
		mv dgdl_${type}\$i.xvg ../;
		cd ../;
	done
	"""

}

process FREE_ENERGY_EST {
/*
	Process executing analysis through pmx to provide a free energy estimate for folded/unfolded structure.
	Input:
		collection for forward transition dgdl xvg files
		collection of reverse transition dgdl xvg files
		type of structure (folded/unfolded)
		name of structure
	Output:
		text file containing free energy
		png of graoh produced from analysis
		.dat file from analysis
*/
	publishDir "${params.pub_dir}/free_energy/${type}/${name}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	tag "${name}"

	input:
	path forward
	path reverse
	val type
	val name

	output:
	path "*.txt"
	path "*.png"
	path "*.dat"

	script:
	"""
	path_f=\$(dirname ${forward[0]})
	path_r=\$(dirname ${reverse[0]})
	pmx analyse -fA \$path_f/dgdl_for* -fB \$path_r/dgdl_rev* -t 298.15 -m bar -o ${type}_${name}_free_energy_est.txt
	"""

}

process CREATE_TRIPEPTIDE {
/*
	Process that makes use of a Python script to create tripeptides.
	Input:
		Tuple containing the index of the mutation, the mutant residue and the native residue
	Output:
		Pdb file of created tripeptide
*/
	publishDir "${params.pub_dir}/created_tripeptides/${native_res}${native_res_index}${mutant}", mode: 'copy', overwrite: false
	container "${simgDir}/gro_pmx_2023.sif"
	tag "native: ${native_res}, mutant: ${mutant}, index: ${index}"

	input:
	tuple val(index), val(mutant), val(native_res)

	output:
	path "*_tripeptide.pdb"

	script:
	"""
	python3 ${PWD}/bin/tripeptide.py -res_native ${native_res} -res_native_index ${index} -res_mutant ${mutant}
	"""


}



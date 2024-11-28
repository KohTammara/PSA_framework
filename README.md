# PSA_framework

The PSA framework is a protein structural analysis pipeline that makes use of NextFlow as its workflow manager. The tools incorporated into this pipeline includes the Rosetta Fixed Backbone design application, the Rosetta SimpleThreadingMover, Maestro - Multi Agent Stability Predition tool, pmx, and Gromacs. The Rosetta tools allow for the implementation of a point mutation on the native structure (Rosetta Fixed Backbone), as well as the implementation of a point mutation on the amino acid sequence of the native protein which is then threaded into the native structure (Rosetta SimpleThreaderMover). The structures produced by the Rosetta processes are then passed to Gromacs, where a molecular dynamics simulation is executed, which helps assess how well the stability of the native structure is maintained with the applied mutation. pmx is a GROMACS compatible tool that can implement mutations as well as analyse GROMACS simulations of the mutants to produce free energy values.This gives insight into the effect of the mutation in terms of the stability of the native structure. Maestro is used to provide a ddG value which allows us to see the destabilising or stabilising effect of the mutation.

## Dependancies

### Nextflow
Nextflow requires Bash 3.2 or later as well as Java 17 (up to version 23).
Instructions on installing Nextflow and its dependancies can be found [here](https://www.nextflow.io/docs/latest/install.html).

**Note** This framework makes use of Nextflow versions 23.10.1 and up.

### Applications (Rosetta, Maestro, pmx, GROMACS)
All applications are containerised and used within the framework. The containerisation recipes can be found in the **Singularity_definition_files** directory. The dependancies of each application is installed within the definition files.
**Note** Rosettas suite download can only occur after obtaining a license. See [Rosetta documentation](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build) for more information.

SingularityCE was used to build the containers, however [Apptainer](https://apptainer.org/docs/user/main/introduction.html) can also be used for the creation and execution of these containers.

### SingularityCE
There are various dependancies required for SingularityCE that can be installed on a Debian based systems as follows:
```
# Ensure repositories are up-to-date
sudo apt-get update
# Install debian packages for dependencies
sudo apt-get install -y \
   autoconf \
   automake \
   cryptsetup \
   fuse2fs \
   git \
   fuse \
   libfuse-dev \
   libglib2.0-dev \
   libseccomp-dev \
   libtool \
   pkg-config \
   runc \
   squashfs-tools \
   squashfs-tools-ng \
   uidmap \
   wget \
   zlib1g-dev
```
GO is also required as a dependancy and the latest version can be downloaded through the [GO official website](https://golang.org/dl/).
A [script](GO_Singularityce_latest_2023_10_18.sh) was created for the purpose of installing both GO and singularity, the versions on this script can be changed to later versions.
More information on installing Singularity can be found on the [SingularityCE official website](https://docs.sylabs.io/guides/main/user-guide/quick_start.html#quick-installation-steps).

## Usage

### Executing the framework

The command to execute the framework is as follows:

`nextflow run main.nf -c nextflow.config -params-file params.json`

`-c` specifies the configuration file

`-params-file` specifies the parameter file

### Input files

Various input files are required for executing this framework depending on the applications used. Each application has its own set of input files.
This framework also makes use of a configuration file as well as a paramete file to provide Nextflow with the required input files.

An example set of input files is available for the mutant M205A for the human Prion protein at ~/test/MUTANTS/M205A_example.

#### Nextflow

The [nextflow configuration](nextflow.config) file will have to be edited to suit the paths of your system, i.e. parts of the absolute paths defined with the nextflow.config file will have to be changed to suit the paths within the host executing the framework. The [params.json](params.json) file would also need to be populated according to the location as well as name of the input files used for the applications within the framework.

#### GROMACS

GROMACS requires a protein structure file as well as the molecular dynamics parameter (MDP) files. The PDB file can be obtained from the Protein Data Bank, and the MDP files can be found in the preproccessed_input firectory. The MDP files for free energy simulations are found in the preprocessed_input/folded_md, which consists of MDP files for the forwards and reverse transitions. The preprocessed_input/md directory houses the MDP files used for a standard MD simulation. The file structure of the MDP directories and files can be seen below:

**Note** The MDP file used for ion generation for the free energy simulations is simply found within the preprocessed_input directory as the file is used for both folded and unfolded free energy simulations. The MDP parameters provided in these files are also specifically for the CHARMM36m forcefield.

```bash
preprocessed_input/
├── folded_md
│   ├── forward
│   │   ├── equil_md
│   │   │   ├── f_enmin.mdp
│   │   │   ├── f_equil.mdp
│   │   │   └── f_npt.mdp
│   │   └── nonequil_md
│   │       └── f_nonequil.mdp
│   └── reverse
│       ├── equil_md
│       │   ├── r_enmin.mdp
│       │   ├── r_equil.mdp
│       │   └── r_npt.mdp
│       └── nonequil_md
│           └── r_nonequil.mdp
├── genion.mdp
├── md
│   ├── em.mdp
│   ├── ions.mdp
│   ├── md.mdp
│   ├── npt.mdp
│   └── nvt.mdp
├── readme.txt
└── unfolded_md
    ├── forward
    │   ├── equil_md
    │   │   ├── f_enmin.mdp
    │   │   ├── f_equil.mdp
    │   │   └── f_npt.mdp
    │   └── nonequil_md
    │       └── f_nonequil.mdp
    └── reverse
        ├── equil_md
        │   ├── r_enmin.mdp
        │   ├── r_equil.mdp
        │   └── r_npt.mdp
        └── nonequil_md
            └── r_nonequil.mdp
```

#### pmx

pmx requires that the PDB file used is preproccessed for the force field that is used. For the purpose of this framework, the CHARMM36m all-atom forcefield was used, and thus the PDB file was preprocessed using CHARMM-gui. More information on this can be found [here](preprocessed_input/readme.txt). A TSV file is also required to specify the mutation to be implemented for pmx. The format is as follows:

```
Mutant  index
<mutant residue>    <index of native residue relative to stucture>
```

The indexing for the PDB structure 1QLX as well as the available forcefield names can be found in [pmx_info.txt](pmx_info.txt)
#### Rosetta

Specific input files are required for each Rosetta application. The same PDB file is used for both applications and is supplied through the parameter `list_of_structs` in [params.json](params.json).
**Note** The PDB file used for rosetta applications does not require preprocessing for forcefield like in the case of pmx.

##### Fixed Backbone application

This application requires a textfile in the format of Rosettas resfile. More information on the syntax and commands for the resfile can be found in the [Rosetta Documentation](https://docs.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles).

To understand how this application works, one can have a look at the [demo](https://docs.rosettacommons.org/demos/latest/public/fixbb_design/README) for Rosettas fixed backbone application and more information about the application itself and how it works can be found [here](https://docs.rosettacommons.org/docs/latest/application_documentation/design/fixbb).

**Note** The usage of this application within this framework has been tested for the simplest usage of mutating a single structure.

##### Simple Threading Mover

This threading application makes use of Rosettas [Movers](https://docs.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/Mover), and threads a mutated sequence into a template structure.
A JSON file is required as input for this application, and it contains the mutant, a start and end position. The positions indicate what portion of the given sequence will be threaded.

#### Meastro

Meastro requires as input a text file which provides the path to the PDB file that will be used as well as the mutation in the following format:

```
<native amino acid><index of native residue>.<chain><{mutant amino acid}>
```

E.g M205.A{A}

### Output

The output of the individual applications are used for free energy simulations to estimate a ddG value in order to assess the effect on stability a mutation has.
Thus the final output from the framework, consists of the output of each application as well as a text file containing the differences in free energies produced by the free energy simulation.

### Limitations

The capabilities of this framework is limited to that of the applications used, to that effect, pmx is limited to mutations that are not charged and does not include Proline, Glycine and Cysteine.

Alongside that, the GROMACS command, `genion`, cannot neutralise both states of the tripeptides obtained from the pmx tripeptide database, which causes errors along the simulation process. 

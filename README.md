# PSA_framework

The PSA framework is a protein structural analysis pipeline that makes use of NextFlow as its workflow manager. The tools incorporated into this pipeline includes the Rosetta Fixed Backbone design application, the Rosetta SimpleThreadingMover, Maestro - Multi Agent Stability Predition tool, pmx, and Gromacs. The Rosetta tools allow for the implementation of a point mutation on the native structure (Rosetta Fixed Backbone), as well as the implementation of a point mutation on the amino acid sequence of the native protein which is then threaded into the native structure (Rosetta SimpleThreaderMover). The structures produced by the Rosetta processes are then passed to Gromacs, where a molecular dynamics simulation is executed, which helps assess how well the stability of the native structure is maintained with the applied mutation. pmx is a GROMACS compatible tool that can implement mutations as well as analyse GROMACS simulations of the mutants to produce free energy values.This gives insight into the effect of the mutation in terms of the stability of the native structure. Maestro is used to provide a ddG value which allows us to see the destabilising or stabilising effect of the mutation.

## Dependancies

### Nextflow
Nextflow requires Bash 3.2 or later as well as Java 17 (up to version 23).
Instructions on installing Nextflow and its dependancies can be found [here](https://www.nextflow.io/docs/latest/install.html).

*Note* This framework makes use of Nextflow versions 23.10.1 and up.

### Applications (Rosetta, Maestro, pmx, GROMACS)
All applications are containerised and used within the framework. The containerisation recipes can be found in the *Singularity_definition_files* directory. The dependancies of each application is installed within the definition files.
*Note* Rosettas suite download can only occur after obtaining a license. See [here](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build) for more information.

SingularityCE was used to build the containers, however [Apptainer](https://apptainer.org/docs/user/main/introduction.html) can also be used for the creation and execution of these containers.

##### SingularityCE
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
GO is also required as a dependancie and the latest version can be downloaded through their [official website](https://golang.org/dl/).
A [script](PSA_framework/GO_Singularityce_latest_2023_10_18) was created for the purpose of installing both GO and singularity, the versions on this script can be changed to later versions.
More information on installing Singularity can be found on their [official website](https://docs.sylabs.io/guides/main/user-guide/quick_start.html#quick-installation-steps).

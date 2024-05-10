PDB format files from the Protein Data Bank or other sources usually require minor modifications with respect to residue 
and atom names before they can be read into CHARMM or Amber, the main modeling packages used by the MMTSB Tool Set. 
CHARMM also requires unique segment IDs at the end of each PDB line.
The pdb files used in this pipeline should be preprocessed according to the forcefield used.

In the case of this sudy, the charmm36m all atom force field was used and preprocessed using CHARMM-GUI for use in pmx and gromacs.
https://cgapi.cc.lehigh.edu/?doc=input/pdbreader

Mutation on input data are implemented within the computational framework and are not implemented with CHARMM-GUI.

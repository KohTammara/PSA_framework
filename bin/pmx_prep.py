from pmx import *
import argparse

parser = argparse.ArgumentParser(description='hybrid file preperation')

parser.add_argument("-pdb", type=str, metavar='/some/path/to/fasta.fasta', help='The path to the protein file', required=True)
parser.add_argument("-res_num", type=int, metavar='120',  help='The positon of the residue to be mutated', required=True)
parser.add_argument("-res_mut", type=str, metavar='eg. F (Phe)',  help='Residue to be mutated to', required=True)
parser.add_argument("-ff", type=str, help='Forcefield to use: amber99sbmut.ff;amber99sb-star-ildn-mut.ff;charmm22starmut.ff;charmm36mut.ff;oplsaamut.ff')
args = parser.parse_args()


# load, mutate, and save the PDB file
m = Model(args.pdb, rename_atoms=True)
m2 = mutate(m=m, mut_resid=args.res_num, mut_resname=args.res_mut, ff=args.ff)
m2.write('mutant.pdb')

# run pdb2gmx
gmx.pdb2gmx(f='mutant.pdb', o='conf.pdb', p='topol.top', ff=args.ff, water='tip3p')

# load topology, fill B states, then write a new topology file
topol = Topology('topol.top')
pmxtop, _ = gen_hybrid_top(topol)
pmxtop.write('newtop.top')
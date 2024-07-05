from pmx.chain import Chain as c
import argparse

parser = argparse.ArgumentParser(description='tripeptide file preperation')
parser.add_argument("-res_native", type=str, metavar='eg. M (Met)', help='Residue abbreviation of orgial residie to be mutated')
parser.add_argument("-res_native_index", type=str, metavar='eg. 81', help='Residue abbreviation of orgial residie to be mutated')
parser.add_argument("-res_mutant", type=str, metavar='eg. F (Phe)',  help='Residue abbreviation of mutant residue', required=True)
args = parser.parse_args()

seq = "G" + args.res_mutant + "G"
chain = c().create(seq)
chain.add_nterm_cap()
chain.add_cterm_cap()
file_name = args.res_native + args.res_native_index + args.res_mutant + "_tripeptide.pdb"
chain.write(file_name)
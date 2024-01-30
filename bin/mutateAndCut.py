#!/usr/bin/env python
import sys
import argparse

parser = argparse.ArgumentParser(description='Mutates and cuts sequence for threading into a template structure with rosetta SimpleThreadingMover')

parser.add_argument("-seq", type=str, metavar='/some/path/to/fasta.fasta', help='The path to the sequence to be threaded stored in a fasta file format', required=True)
parser.add_argument("-mutation", type=str, metavar='"M205A"',  help='Mutation to be applied to seqeunce, position should correspond to original numbering as given in your input, format as "M205A"(M-Methionine at position 205 to A-Alanine)', required=True)
parser.add_argument("-start_position", type=int, default=1, help='Optional argument to define a start position if it is not 1')
parser.add_argument("-end_position", default=False, help='Optional argument to define a end position if required, can either be an integer or False')
args = parser.parse_args()

end_file = open("new_seq.fasta", "w")

def mutate_and_cut(file_path, mutant_with_position, start_pos=1, end_pos=False):
	# Read in the sequence from a fasta file
	f = open(file_path,"r")
	lines = f.readlines()
	seq = []
	for line in lines:
		if '>' not in line:
			seq.append(line.strip("\n"))
	fasta_seq = "".join(seq)

	#Mutate the seqeunce at the given position
	existing_residue = mutant_with_position[0:1]
	mutant_residue = mutant_with_position[len(mutant_with_position)-1:len(mutant_with_position)]
	position = int(mutant_with_position[1:len(mutant_with_position)-1])
	if end_pos == "False":
		end_pos = False
	else:
		end_pos = int(end_pos)
	if existing_residue == fasta_seq[position-1]:
		# print(f'{existing_residue} and then {position} {fasta_seq[position-1]}')
		seq = list(fasta_seq)
		seq[position-1] = mutant_residue
		mutated_seq = ''.join(seq)
		cut_seq = ''
		#Cut the sequence according to the given positions
		if end_pos is not False:
			cut_seq = mutated_seq[start_pos-1:end_pos]
		else:
			cut_seq = mutated_seq[start_pos-1:]
		print(cut_seq)
		end_file.write(str(cut_seq))
		end_file.close()
	else:
		print('The position provided does not match the original amino acid that is to be mutated.', file=sys. stderr)
		sys.exit(1)


mutate_and_cut(args.seq, args.mutation, args.start_position, args.end_position)


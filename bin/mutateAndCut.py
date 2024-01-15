import sys
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
	if existing_residue == fasta_seq[position-1]:
		# print(f'{existing_residue} and then {position} {fasta_seq[position-1]}')
		seq = list(fasta_seq)
		seq[position-1] = mutant_residue
		mutated_seq = ''.join(seq)

		#Cut the sequence according to the given positions
		if end_pos is not False:
			cut_seq = mutated_seq[start_pos-1:end_pos]
		else:
			cut_seq = mutated_seq[start_pos-1:]
		
		return fasta_seq,cut_seq 
	else:
		print('The position provided does not match the original amino acid that is to be mutated.', file=sys. stderr)
		sys.exit(1)


if __name__ == "__main__":
	fasta_seq,cut_seq = mutate_and_cut('Q53YK7.fa', 'M207A', 1)
	print(f'behold the seqeunce to not use in rosetta fixbb {fasta_seq}\nbehold the seqeunce to     use in rosetta fixbb {cut_seq}')


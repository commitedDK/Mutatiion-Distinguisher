#!/usr/bin/python3



#this code takes 3 files as input:- gff file, reference genome file, and vcf file which contains mutations, the task is to figure out that point mutations are synonymous or non-synonymous. 
# sys.argv[1] = Gallus.gff3
# sys.argv[2] = Gallus.fa
# sys.argv[3] = variant.vcf
# sys.argv[4] = output.txt


import sys
def main():# this function call all the functions used
	Gff = sys.argv[1]
	Reference = sys.argv[2]
	Vcf = sys.argv[3]

	Upperlimit_list, Lowerlimit_list, Phase_list = make_list_of_CDS(Gff)
	seq = extract_seq(Reference)
	seq_list = list(seq)
	CDS_dict, CDS_seq_list = make_dict_of_CDSseq(Upperlimit_list,Lowerlimit_list,Phase_list,seq)
	Position_list, Mutation_list = make_list_of_Positions(Vcf)
	
	Position_list =list(map(int, Position_list))
	mutated_seq_list = replace_positions_in_ref_genome(seq_list,Position_list, Mutation_list)
	mutated_seq = ''.join(mutated_seq_list)

	mutated_CDS_dict, mutated_CDS_seq_list = make_dict_of_CDSseq(Upperlimit_list,Lowerlimit_list,Phase_list,mutated_seq)
	
	dict_check = Check_type_of_mutation(CDS_dict, mutated_CDS_dict)
	feature_dict = make_list_of_Positions_in_genome(dict_check)
	print(feature_dict)


def make_list_of_CDS(input_file): # this function considers Coding regions only and make lists of start positions, end positions and phase of the CDS.
	with open(input_file, 'r') as fin:
		Upperlimit_list = []
		Lowerlimit_list = []
		Phase_list = []
		for line in fin:
			if 'CDS' in line :
				line = line.strip()
				Start = line.split('\t')[3]
				End = line.split('\t')[4]
				Phase = line.split('\t')[7]
				Upperlimit_list.append(int(End))
				Lowerlimit_list.append(int(Start))
				Phase_list.append(int(Phase))

		return Upperlimit_list, Lowerlimit_list, Phase_list


def extract_seq(input_file): # this function converts multiple lines of sequence into one line.
	with open(input_file,'r') as file_fna:
		seq_list = []
		ids = []
		seq = ''
		for line in file_fna:
			line = line.rstrip()
			if line.startswith('>'):
				if seq != '':
					seq_list.append(seq)
					seq = ''
				ids.append(line)
			else:
				seq += line
		seq_list.append(seq)


	return seq

def make_dict_of_CDSseq(Upperlimit_list, Lowerlimit_list,Phase_list,seq) : # this function creates dictionaries with key containing information of phase, start postion and end position of CDS and value as sequence of CDS.
	CDS_seq_list = []
	CDS_dict = {}
	for i in range(len(Upperlimit_list)):
		upper = Upperlimit_list[i]
		lower = Lowerlimit_list[i]
		phase = Phase_list[i]
		CDS_seq = seq[lower-1:upper-1]
		CDS_seq_list.append(CDS_seq)
		key = 'CDS' + str(phase) + ' ' + str(lower) + ':' + str(upper)
		CDS_dict[key] = CDS_seq

	return CDS_dict, CDS_seq_list

def make_list_of_Positions(input_file): # this function is used to make dictionary of positions, with position as key and mutated nucleotide as value.
	with open(input_file, 'r') as fin:
		Position_list = []
		Position_dict = {}
		Mutation_list = []
		for line in fin:
			if 'dbSNP_150;TSA=SNV' in line :
				line = line.strip()
				Position = line.split('\t')[1]
				Position = int(Position) - 1
				Ref = line.split('\t')[3]
				Alt = line.split('\t')[4]
				if ',' in Alt :
					Alt = Alt.split(',')[0]
				elif len(Alt) > 1 :
					Alt = Alt.split('')[0]
				Position_list.append(Position)
				Mutation_list.append(Alt)
				Position_dict[Position] = Alt
		return Position_list, Mutation_list

def replace_positions_in_ref_genome(to_modify,indexes,replacements) :
    for a_element, m_element in zip(indexes, replacements):
        to_modify[a_element] = m_element

    return(to_modify)

def Check_type_of_mutation(CDS_dict, mutated_CDS_dict):
	dict_check = {}
	amino_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S",
            "TCA":"S","TCG":"S","TAT":"Y", "TAC":"Y", "TAA":"STOP",
             "TAG":"STOP","TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
             "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L","CCT":"P", 
             "CCC":"P", "CCA":"P", "CCG":"P","CAT":"H", "CAC":"H", 
             "CAA":"Q", "CAG":"Q","CGT":"R", "CGC":"R", "CGA":"R", 
             "CGG":"R","ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
             "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T","AAT":"N", 
             "AAC":"N", "AAA":"K", "AAG":"K","AGT":"S", "AGC":"S", 
             "AGA":"R", "AGG":"R","GTT":"V", "GTC":"V", "GTA":"V", 
             "GTG":"V","GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
             "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E","GGT":"G", 
             "GGC":"G", "GGA":"G", "GGG":"G", "NTC":"unknown"}
	for key in CDS_dict:
		list_check = []
		value = CDS_dict[key] # saving the values of dictionary into variables to easily work on.
		mutated_value = mutated_CDS_dict[key]
		
		if (value != mutated_value): # checks if values of CDS_dict and mutated_CDS_dict is different.

			if key.startswith('CDS0') :
				startindex = 0
				#print(key,startindex)
			elif key.startswith('CDS1') :
				startindex = 1
	            #print(key,startindex)
			elif key.startswith('CDS2') :  
				startindex = 2
	            #print(key,startindex)
			for i in range(startindex,len(mutated_value)) : # 0,1,2,3,4,5,... basepair index of mutated exon
				x=3*((i-startindex)//3)+startindex  # 0,0,0, 3,3,3, 6,6,6, 9,9,9,... index of codon start positions in exon
				#print('x',x,'i',i)
				#codon = '' 
				#mutated_codon = ''
				codon = value[x:x+3] 
				mutated_codon = mutated_value[x:x+3]
	            
				if (len(mutated_codon) != 3) :
					break

				if (value[i]!= mutated_value[i]) : # If basepair is mutated, write it to the list

					if (amino_dict[codon] != amino_dict[mutated_codon]) : # makes different amino acid
						list_check.append((i,'non-syn',codon,mutated_codon))
						#print(list_check)
					else : # makes same amino acid
						list_check.append((i,'syn',codon,mutated_codon))
					dict_check[key]=list_check

	return(dict_check)

def make_list_of_Positions_in_genome(Dict): # this function is used to make dictionary of positions, with position as key and mutated nucleotide as value.
	feature_dict = {}
	for Key,value in Dict.items() :			
		Key = Key.split('\t')[0]
		Key = Key.split(' ')[1]
		Key = int(Key.split(':')[0])
		for v in value :

			Value = v[0]
			Type = v[1]
			position = (Key + Value - 1)
			feature_dict[position] = Type
			#print(position, Type)
	for index, value in sorted(feature_dict.items()):
		#print(index,value)

		return feature_dict 
	

def write_to_file(output_file,dict_check) :
    with open(output_file, 'w') as fout:
        print(dict_check, file = fout)

main()


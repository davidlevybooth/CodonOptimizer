#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from CAI import RSCU
from dnachisel import (
	DnaOptimizationProblem,
	random_protein_sequence,
	reverse_translate,
	CodonOptimize,
	EnforceTranslation,
	AvoidPattern,
	AvoidChanges,
	EnforceSequence,
	EnforceGCContent,
)

ap = argparse.ArgumentParser()
ap.add_argument("-e", "--highexp", required=False,
	help="path to the input fa file", type=str)
ap.add_argument("-i", "--inputfasta", required=True,
	help="path to the gene to be optimized (fasta file)", type=str)
ap.add_argument("-t", "--taxid", required=False,
	help="Host organism NCBI Taxonomic ID", type=str)
ap.add_argument("-c", "--genetic_code", required=False,
	help="Genetic code # for organism", type=int, default=11)
ap.add_argument("-r", "--report", required=False,
	help="Path for output report", type=str)
ap.add_argument("-o", "--output", required=True,
	help="Path for fasta ouput", type=str)
ap.add_argument("-m", "--method", required=False,
	help="Method for codon optimization", type=str, default="match_codon_usage")
ap.add_argument("-p", "--protein",
	help="Flag if protein input", action='store_true')

args = vars(ap.parse_args())

input_path = args["highexp"]
fasta_path = args["inputfasta"]
taxid = args["taxid"]
genetic_code = args["genetic_code"]
report_path = args["report"]
output_path = args["output"]
protein_flag = args["protein"]

## Codon table
#``{'*': {"TGA": 0.112, "TAA": 0.68}, 'K': ...}``
codon_table_11={"*": {"TGA":1.0, "TAA":1.0, "TAG":1.0},
"T": {"ACA":0, "ACC":0, "ACG":0, "ACT":0},
"R": {"AGA":0, "AGG":0, "CGA":0, "CGC":0, "CGG":0, "CGT":0},
"A": {"GCA":0, "GCC":0, "GCG":0, "GCT":0},
"C": {"TGC":0, "TGT":0},
"D": {"GAC":0, "GAT":0},
"E": {"GAA":0, "GAG":0},
"F": {"TTC":0, "TTT":0},
"G": {"GGA":0, "GGC":0, "GGG":0, "GGT":0},
"H": {"CAC":0, "CAT":0},
"I": {"ATA":0, "ATC":0, "ATT":0},
"K": {"AAA":0, "AAG":0},
"L": {"CTA":0, "CTC":0, "CTG":0, "CTT":0, "TTA":0, "TTG":0},
"M": {"ATG":1.0},
"N": {"AAC":0, "AAT":0},
"P": {"CCA":0, "CCC":0, "CCG":0, "CCT":0},
"Q": {"CAA":0, "CAG":0},
"S": {"AGC":0, "AGT":0, "TCA":0, "TCC":0, "TCG":0, "TCT":0},
"V": {"GTA":0, "GTC":0, "GTG":0, "GTT":0},
"W": {"TGG":0},
"Y": {"TAC":0, "TAT":0}}

if genetic_code == 11:
	codon_table = codon_table_11
else:
	print("\ngenetic codes other than 11 (Bacterial, Archaeal) not supported")


#Import target gene
#To do: refactor to process multiple input sequences
gene_object = SeqIO.parse(fasta_path, "fasta")
for dna_seq in gene_object:
			dna_id = dna_seq.id
			print("\nImporting target gene " + dna_id)
			dna = str(dna_seq.seq)
			gene = dna


if input_path and not taxid:
	print("\ngene list and taxonomic ID both provided. Defaulting to gene list")

#Import gene list for RSCU calculation
if input_path and not taxid:
	print("\nImporting genes for RSCU calculation")
	seq_list = []
	counter = 0
	n_count = 0
	seq_object = SeqIO.parse(input_path, "fasta")
	for seqs in seq_object:
				seq_id = seqs.id
				seq = str(seqs.seq)
				seq_list.append(seq)
				counter += 1
				nn = len(seq)
				n_count += nn

	print("\n" + str(counter) + " genes imported containing " + str(n_count) + " nucleotides")
	print("\nCalculating RSCU for imported genes\n")
	try:
		RSCU_list = RSCU(seq_list)
	except:
		print("\nEXCEPTION: RSCU could not be caluclated for imported genes")
	print("\nParsing RSCU for codon optimization")
	#To do: Create RSCU parser function
	for k, v in codon_table_11.items():
		for k2, v2 in v.items():
			if k2 in RSCU_list:
				codon_table_11[k][k2] = RSCU_list[k2]

	print("\nOptimizing codons for input gene list")
	#Read gene fasta sequence and initiate optimizer
	#To do: populate constraints in function to prevent duplication
	if not protein_flag:
		problem = DnaOptimizationProblem(
			sequence=gene,
			constraints=[
				AvoidPattern("BsmBI_site", "BamHI"),
				EnforceTranslation(),
				AvoidChanges(location=(0, 2)),
				#EnforceSequence(sequence = "ATG", location=(0, 2)),
				EnforceGCContent(mini=0.35, maxi=0.65, window=50), #TWIST: 25% and 65% GC
			],
			objectives=[CodonOptimize(codon_usage_table=codon_table_11)],
		)
	if protein_flag:
		gene = reverse_translate(gene)
		problem = DnaOptimizationProblem(
		sequence=gene,
		constraints=[
			AvoidPattern("BsmBI_site", "BamHI"),
			EnforceTranslation(),
			EnforceGCContent(mini=0.35, maxi=0.65, window=50), #TWIST: 25% and 65% GC
		],
		objectives=[CodonOptimize(codon_usage_table=codon_table_11)],
	)


if taxid and not input_path:
	print("\nOptimizing codons for taxonomic ID: " + taxid)
	#Read gene fasta sequence and initiate optimizer
	if not protein_flag:
		problem = DnaOptimizationProblem(
			sequence=gene,
			constraints=[
				#EnforceSequence(sequence = "ATG", location=(0, 2)),
				AvoidChanges(location=(0, 2)),
				AvoidPattern("BsmBI_site", "BamHI"),
				EnforceTranslation(),
				EnforceGCContent(mini=0.35, maxi=0.65, window=50), #TWIST: 25% and 65% GC

			],
			objectives=[CodonOptimize(species=taxid)],
		)
	if protein_flag:
		gene = reverse_translate(gene)
		problem = DnaOptimizationProblem(
		sequence=gene,
		constraints=[
			EnforceTranslation(),
			AvoidPattern("BsmBI_site", "BamHI"),
			EnforceGCContent(mini=0.35, maxi=0.65, window=50), #TWIST: 25% and 65% GC
		],
		objectives=[CodonOptimize(species=taxid)],
	)

#Output and reporting
print("\nBefore optimization:")
print(problem.constraints_text_summary())
print(problem.objectives_text_summary())

problem.resolve_constraints(final_check=True)

if report_path:
	print("\nSaving report " + report_path)
	problem.optimize_with_report(target=report_path)
else:
	problem.optimize()

print("\nAfter optimization:")
print(problem.constraints_text_summary())
print(problem.objectives_text_summary())

print("\nWriting optimized sequence as fasta: " + output_path)


dna_id_new=dna_id+"_opt"
dna_seq_new=Seq(problem.sequence, generic_dna)

records = (SeqRecord(dna_seq_new, id=dna_id_new, description=""))
with open(output_path, "w") as output_handle:
	SeqIO.write(records, output_handle, "fasta")

protein_seq_new=dna_seq_new.translate(table=11, cds=True)
protein_records=(SeqRecord(protein_seq_new, id=dna_id_new, description=""))

protein_path=output_path.split(".")[0]+".faa"
with open(protein_path, "w") as output_handle:
	SeqIO.write(protein_records, output_handle, "fasta")

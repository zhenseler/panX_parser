'''
Usage: Format panX data into a human readable format, which can be used as input for Roary

python panX_parser.py <input_folder> <genome_list_file> <output_filename>

input_folder: folder containing panX output
genome_list_file: a return-seperated file with all genome names
output_filename: name you want the output file to have

-Z
'''

# Imports
from sys import argv
import gzip
import os
from collections import Counter

# Args
input_folder = argv[1]
if not input_folder.endswith("/"):
	input_folder += "/"
genome_list_file = argv[2]
output_filename = argv[3]

out = open(input_folder + output_filename, "w")

output_header = ("\t").join(["File", "Gene", "Non-unique Gene name", "Annotation", "No. isolates", \
"No. sequences", "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", \
"Accessory Fragment", "Accessory Order with Fragment", "QC", "Min group size", \
"Max group size", "Avg group size"])

# Get all input filenames
encoded_folder = os.fsencode(input_folder)
filenames = []
for filename in os.listdir(encoded_folder):
	filename = os.fsdecode(filename)
	if filename.endswith(".fa.gz"):
		filenames.append(filename)
filenames.sort()

# Make genome list (list of all input genomes
genome_list = []
for l in open(genome_list_file, "U"):
	genome_list.append(l.strip("\n"))
genome_list.sort()

output_header_tail = ("\t").join(genome_list)
out.write(output_header + "\t" + output_header_tail + "\n")

# Loop through files
for file in filenames:
	print("\nRunning one file " + file)
	full_filepath = input_folder + file	

	# Count number of lines in file
	file_line_counter = 0
	for l in gzip.open(full_filepath, "rt"):
		file_line_counter += 1
	# Keeping track of sequence lengths for min and max
	seq_lengths = []
	# list of all instances of each strain in that file
	strain_list = []
	# Used to associate locus_tag with strain
	locus_list = []
	seq = ""
	presence_absence_list = []
	
	for i,l in enumerate(gzip.open(full_filepath, "rt")):
		if l.startswith(">"):			
			# Unpack info stored from previous sequence
			if len(seq) != 0:
				seq_lengths.append(len(seq)-seq.count("-"))
			seq = ""
			
			hyphen_split = l.split("-")
			
			strain = hyphen_split[0].lstrip(">")
			strain_list.append(strain)
			
			locus_tag = hyphen_split[1]
			locus_list.append(locus_tag)

			annotation = ("-").join(hyphen_split[3:])
			annotation = annotation.split("_")
			if annotation[1].isdigit():
				geneID = annotation[0]
				annotation = "_".join(annotation[2:])
			else:
				geneID = annotation[0]
				annotation = "_".join(annotation[1:])
				
			# Fields that remain blank in this format
			non_unique_gene_name = ""
			genome_fragment = ""
			order_within_fragment = ""
			accessory_fragment = ""
			accessory_order_within_fragment = ""
			QC = ""
			
		else:
			seq += l.strip("\n")
		# Get length of last sequence
		if i == file_line_counter - 1:
			seq_lengths.append(len(seq)-seq.count("-"))
			
	number_of_isolates = len(set(strain_list))
	number_of_seqs = len(strain_list)
	avg_number_of_seqs_per_isolate = len(strain_list)/float(len(set(strain_list)))
	
	# If the same genome appears twice, combines all instances in file with hyphens
	for genome in genome_list:
		locus_tag = ""
		for i,strain in enumerate(strain_list):
			if strain == genome:
				if locus_tag == "":
					locus_tag = locus_list[i]
				else:
					locus_tag += "," + locus_list[i] # Separator for when one genome has multiple loci
			if i == len(strain_list) - 1:
				presence_absence_list.append(locus_tag)

	minimum_seq_length = min(seq_lengths)
	maximum_seq_length = max(seq_lengths)
	avg_seq_length = sum(seq_lengths)/float(len(seq_lengths))
	
	filename_for_output = file.rstrip(".fa.gz")
	
	output_front = ("\t").join([filename_for_output, geneID, non_unique_gene_name,\
	 annotation.strip("<unknown description>\n"), str(number_of_isolates), \
	 str(number_of_seqs), str(avg_number_of_seqs_per_isolate), genome_fragment, \
	 order_within_fragment, accessory_fragment, accessory_order_within_fragment,\
	  QC, str(minimum_seq_length), str(maximum_seq_length), str(avg_seq_length)])
	 
	output_back = ("\t").join(presence_absence_list)
	output_final = ("\t").join([output_front, output_back])+"\n"
	
	out.write(output_final)
	
out.close()

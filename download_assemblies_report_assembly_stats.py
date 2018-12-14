# Program downloads GenBank assemblies for the genera of protist pathogens of interest.
# Program collects download statistics for each protist pathogen.
# Program concatenates assemblies into single fasta database.
# Program reports length (in bp) of each assembly.

# User must provide output directory
# Example command: python download_assemblies_report_assembly_stats.py /home/user/output_directory/

import os
import sys
from Bio import SeqIO

output_directory = sys.argv[1]

# Program downloads the assembly metadata provided by GenBank 
os.system('wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt')
metadata = "assembly_summary_genbank.txt"

# Genera with known protist pathogens 
parasites = {"Karlodinium":{},
"Prorocentrum":{},
"Chattonella":{},
"Dinophysis":{},
"Alexandrium":{},
"Pyrodinium":{},
"Pseudo-nitzschia":{},
"Azidinium":{},
"Karenia":{},
"Toxoplasma":{},
"Giardia":{},
"Entamoeba":{},
"Cryptosporidium":{},
"Cyclospora":{},
"Balamuthia":{},
"Naegleria":{},
"Acanthamoeba":{},
"Trypanosoma":{},
"Gambierdiscus":{},
"Plasmodium":{},
"Leishmania":{},
"Blastocystis":{},
"Babesia":{},
"Balamuthia":{},
"Balantidium":{},
"Cystoisospora":{},
"Dientamoeba":{},
"Sappinia":{},
"Sarcocystis":{},
"Trichomonas":{},
"Trypanosoma":{}}

stats = {}

for h,i in enumerate(open(metadata)):
	if h != 0:
		# parse metadata file
		tmp = i.strip().split('\t')
		taxaName = tmp[7]
		ftp = tmp[19]
		genome = ftp+"/"+ftp.split("/")[-1]+"_genomic.fna.gz"
		# Download if genera of protist pathogen in taxonomic name of assembly
		for k,v in parasites.items():
			if k.upper() in taxaName.upper():
				os.system("wget -nc -O "+output_directory+" "+genome)
				# Record number of assemblies downloaded for each protist pathogen genera and species
				if taxaName in parasites[k]:
					parasites[k][taxaName] += 1
				else:
					parasites[k][taxaName] = 1
				break


# Write out download statistics report
with open(output_directory+'download_statistics.txt','w') as out:
	for k,v in sorted(parasites.items()):
		if v == {}:
			out.write(k+"\t0\n")
		for z,w in v.items():
			out.write(k+"\t"+str(z)+"\t"+str(w)+"\n")
	out.write("\n\n")
	for k,v in sorted(parasites.items()):
		out.write(k+"\t"+str(sum(v.values())+"\n")


# Write to report the bp length of each assembly
with open(output_directory+'protist_pathogens_genomic.fasta','w') as out, open(output_directory+'assembly_stats.txt','w') as out2:
	assembly_length = {}
	for i in os.listdir(output_directory):
		for j in SeqIO.parse(output_directory+i,'fasta'):
			if i in assembly_length:
				assembly_length[i] += len(i.seq)
			else:
				assembly_length[i] = len(i.seq)
			out.write(">"+str(i)+str(j.description).replace(' ','_')+"\n"+str(i.seq)+"\n")
	for k,v in assembly_length.items():
		out2.write(k+"\t"+str(v)+"\n")
		


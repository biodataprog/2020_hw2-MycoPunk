#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re
import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")
    
#step 2: print the number of 'genes'
counter = 0
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff_in = csv.reader(fh,delimiter="\t")
    for row in gff_in:
        if "gene" in row:
            counter +=1
#    print(counter)
    print("there are ",counter,"genes in the genome")

#step 3: Compute the total length of the genes (length is the END - START)
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff_in = csv.reader(fh,delimiter="\t") 
    list_of_length = []
    for row in gff_in:
        if "gene" in row:
#            print(int(row[4]) - int(row[3]))
            list_of_length.append( int(row[4]) - int(row[3]))
    print("the total length of all genes is ",sum(list_of_length)," BP")

#step 4: use fasta file to print the total length of the genome
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"
with gzip.open(fasta,"rt") as fasta_in:	
    for seq_record in SeqIO.parse(fasta_in, "fasta"):
        print("the total length of the genome is ",len(seq_record), " BP")

#step 5: print the % of the genome that is coding
#print(float(sum(list_of_length))/float(len(seq_record)))*100," % of the genome is coding")
print(float(sum(list_of_length))/float(len(seq_record))*100," % of the genome is coding")

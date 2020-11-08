#!/usr/bin/env python3

#load packages
import os,gzip,itertools,csv,re
import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

#input data
url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

#1. The total number of genes in each species.
count = 0
with gzip.open(file1,"rt") as fasta_in:
    for record in SeqIO.parse(fasta_in, "fasta"):
        count += 1
print("Salmonella enterica has ",count," genes")


count2 = 0
with gzip.open(file2,"rt") as fasta_in:
    for record in SeqIO.parse(fasta_in, "fasta"):
        count2 += 1
print("Mycobacterium tuberculosis has ",count2," genes")


#2. Total length of these gene sequences for each file
Se_gene_lengths = []
with gzip.open(file1,"rt") as fasta_inSe:
    for sequence in SeqIO.parse(fasta_inSe, "fasta"):
        Se_gene_lengths.append(len(sequence.seq))    
print("The total length of all genes in Salmonella enterica is",sum(Se_gene_lengths)," BP")

Mt_gene_lengths = []
with gzip.open(file2,"rt") as fasta_inMt:
    for sequence in SeqIO.parse(fasta_inMt, "fasta"):
        Mt_gene_lengths.append(len(sequence.seq))
print("The total length of all genes in Mycobacterium tuberculosis is",sum(Mt_gene_lengths)," BP")


#3. The G+C percentage for the whole dataset (eg the frequency of G + the frequency of C)
Se_G_orC = ''
with gzip.open(file1,"rt") as fasta_inSe:
    for sequence in SeqIO.parse(fasta_inSe, "fasta"):
        Se_G_orC += sequence.seq        
print("The total GC content of Salmonella enterica is",GC(Se_G_orC))


Mt_G_orC = ''
with gzip.open(file2,"rt") as fasta_inMt:
    for sequence in SeqIO.parse(fasta_inMt, "fasta"):
        Mt_G_orC += sequence.seq
print("The total GC content of Mycobacterium tuberculosis is",GC(Mt_G_orC))


print("The combined total GC of the dataset is ",GC(Se_G_orC + Mt_G_orC)) 


#4. Total number codons in each genome.
def split_str(str, codon_size):
   return [str[i:i+codon_size] for i in range(0, len(str), codon_size)]

Se_split_codons = ''
with gzip.open(file1,"rt") as fasta_inSe:
    for sequence in SeqIO.parse(fasta_inSe, "fasta"):
        Se_split_codons += str(sequence.seq)
listed_codons_Se = split_str(Se_split_codons, 3)
print("The total number of codons in Salmonella enterica is ",len(listed_codons_Se))

Mt_split_codons = ''
with gzip.open(file2,"rt") as fasta_inMt:
    for sequence in SeqIO.parse(fasta_inMt, "fasta"):
        Mt_split_codons += str(sequence.seq)
listed_codons_Mt = split_str(Mt_split_codons, 3)
print("The total number of codons in Mycobacterium tuberculosis is ",len(listed_codons_Mt))


#get uninque values of listed_codons, a number of number in each category
unique_codons_Se =list(set(listed_codons_Se))
unique_codons_Mt =list(set(listed_codons_Mt))
#print(unique_codons_Se)
#print(len(unique_codons_Se))
#print(unique_codons_Mt)
#print(len(unique_codons_Se))

codon_dictionary_Se ={}
for k in unique_codons_Se:
    codon_dictionary_Se[k] = listed_codons_Se.count(k)

codon_dictionary_Mt ={}
for k in unique_codons_Mt:
    codon_dictionary_Mt[k] = listed_codons_Mt.count(k)

#5. Print out table with three columns: Codon, Frequency in Sp1, Frequency in Sp2
print("Codon", "\t", "Frequency in Salmonella enterica", "\t", "Frequency in Mycobacterium tuberculosis")
for k in codon_dictionary_Se:
        print(k,"\t", codon_dictionary_Se[k], "\t",codon_dictionary_Mt[k] )

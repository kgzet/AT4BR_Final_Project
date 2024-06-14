# Kinga Zajdel
# AT4BR final project

from Bio import Entrez, SeqIO
from collections import Counter

# preparing variables for future use
counting_CDS = 0
segments_CDS = []
codons_list = []
counted_codons = {}
list_of_counted_aminoacids = {}
genome_id = ["OR682571"]

# table for translating the codons into amino acids
table_of_codons = {
    'AAA': 'K',    # Lysine
    'AAC': 'N',    # Asparagine
    'AAG': 'K',    # Lysine
    'AAT': 'N',    # Asparagine
    'ACA': 'T',    # Threonine
    'ACC': 'T',    # Threonine
    'ACG': 'T',    # Threonine
    'ACT': 'T',    # Threonine
    'AGA': 'R',    # Arginine
    'AGC': 'S',    # Serine
    'AGG': 'R',    # Arginine
    'AGT': 'S',    # Serine
    'ATA': 'I',    # Isoleucine
    'ATC': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ATT': 'I',    # Isoleucine
    'CAA': 'Q',    # Glutamine
    'CAC': 'H',    # Histidine
    'CAG': 'Q',    # Glutamine
    'CAT': 'H',    # Histidine
    'CCA': 'P',    # Proline
    'CCC': 'P',    # Proline
    'CCG': 'P',    # Proline
    'CCT': 'P',    # Proline
    'CGA': 'R',    # Arginine
    'CGC': 'R',    # Arginine
    'CGG': 'R',    # Arginine
    'CGT': 'R',    # Arginine
    'CTA': 'L',    # Leucine
    'CTC': 'L',    # Leucine
    'CTG': 'L',    # Leucine
    'CTT': 'L',    # Leucine
    'GAA': 'E',    # Glutamate
    'GAC': 'D',    # Aspartate
    'GAG': 'E',    # Glutamate
    'GAT': 'D',    # Aspartate
    'GCA': 'A',    # Alanine
    'GCC': 'A',    # Alanine
    'GCG': 'A',    # Alanine
    'GCT': 'A',    # Alanine
    'GGA': 'G',    # Glycine
    'GGC': 'G',    # Glycine
    'GGG': 'G',    # Glycine
    'GGT': 'G',    # Glycine
    'GTA': 'V',    # Valine
    'GTC': 'V',    # Valine
    'GTG': 'V',    # Valine
    'GTT': 'V',    # Valine
    'TAA': '*',    # Stop
    'TAC': 'Y',    # Tyrosine
    'TAG': '*',    # Stop
    'TAT': 'Y',    # Tyrosine
    'TCA': 'S',    # Serine
    'TCC': 'S',    # Serine
    'TCG': 'S',    # Serine
    'TCT': 'S',    # Serine
    'TGA': '*',    # Stop
    'TGC': 'C',    # Cystein
    'TGG': 'W',    # Tryptophane
    'TGT': 'C',    # Cysteine
    'TTA': 'L',    # Leucine
    'TTC': 'F',    # Phenylalanine
    'TTG': 'L',    # Leucine
    'TTT': 'F'    # Phenylalanine
}

# Entrez from Biopython enables downloading data from such sites like NCBI
# here I need nucleotide sequence accompanied by information where are the coding sequences
Entrez.email = ""
handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
# records = SeqIO.parse(handle, "gb")
sequence_info = SeqIO.read(handle, "gb")

# the whole sequence in which I will look for CDS from start to end
whole_sequence = sequence_info.seq

for feat in sequence_info.features:
    # I'm looking for lines with the coding sequence
    if feat.type == 'CDS':
        # counting how much of methionine codons are actually start codons - they will not be counted into "bias"
        # start codons are obligatory, there is no place for the usage bias
        counting_CDS += 1
        # the list of CDS' starts and ends; one list element = one CDS
        segments_CDS.append([int(feat.location.start), int(feat.location.end)])
# print(segments_CDS)
# print(whole_sequence[687], whole_sequence[6057], '\n')
i = 0
# !!
# segments_CDS = []

# loop for checking every pair of start and stop
# "s" stands for single CDS
for s in segments_CDS:
    # print(s[0],"\n",s[1])
    # print(s)
    # "i" is an index according to whole sequence, in range from the start to the end of the CDS
    # range(start, stop, step)
    for i in range(s[0], s[1], 3):
        # building the codon = three subsequent nucleotides (string)
        single_codon = whole_sequence[i] + whole_sequence[i+1] + whole_sequence[i+2]
        # building the list of codons (string)
        codons_list.append(single_codon)

# print(codons_list)
# Counter builds a dictionary with number of occurrences as values
counted_codons = dict(Counter(codons_list))
# print(counted_codons, '\n')
# now I need to reduce the number of Methionines by how many CDS I have = how many start codons
counted_codons['ATG'] -= counting_CDS
# print(counted_codons)

# adding the amino acid shortcut
# new dictionary with key = codon and value = list of two elements: amino acid shortcut and number of occurrences
for k in counted_codons.keys():
    for a in table_of_codons.keys():
        if k == a:
            list_of_counted_aminoacids[table_of_codons[a]] = [k, counted_codons[k]]
            # list_of_counted_aminoacids[k] = [table_of_codons[a], counted_codons[k]]

# print(list_of_counted_aminoacids)

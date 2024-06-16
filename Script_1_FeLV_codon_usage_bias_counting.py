# Kinga Zajdel
# AT4BR final project

from Bio import Entrez, SeqIO
from collections import Counter
import csv

# preparing variables for future use
genome_id = ["OR682571"]        # put here IDs of genomes to examine
counting_CDS = 0        # how many CDS there are, this will be used to count Methionines
segments_CDS = []       # coding sequences sliced from the whole genome
codons_list = []        # list of codons existing in the NCBI sequence
list_of_counted_aminoacids = []     # final list to export into .csv file

# dictionary of codon : how many; setting a list of every possible triplet with starting value = 0
counted_codons = {
    'AAA': 0,
    'AAC': 0,
    'AAG': 0,
    'AAT': 0,
    'ACA': 0,
    'ACC': 0,
    'ACG': 0,
    'ACT': 0,
    'AGA': 0,
    'AGC': 0,
    'AGG': 0,
    'AGT': 0,
    'ATA': 0,
    'ATC': 0,
    'ATG': 0,
    'ATT': 0,
    'CAA': 0,
    'CAC': 0,
    'CAG': 0,
    'CAT': 0,
    'CCA': 0,
    'CCC': 0,
    'CCG': 0,
    'CCT': 0,
    'CGA': 0,
    'CGC': 0,
    'CGG': 0,
    'CGT': 0,
    'CTA': 0,
    'CTC': 0,
    'CTG': 0,
    'CTT': 0,
    'GAA': 0,
    'GAC': 0,
    'GAG': 0,
    'GAT': 0,
    'GCA': 0,
    'GCC': 0,
    'GCG': 0,
    'GCT': 0,
    'GGA': 0,
    'GGC': 0,
    'GGG': 0,
    'GGT': 0,
    'GTA': 0,
    'GTC': 0,
    'GTG': 0,
    'GTT': 0,
    'TAA': 0,
    'TAC': 0,
    'TAG': 0,
    'TAT': 0,
    'TCA': 0,
    'TCC': 0,
    'TCG': 0,
    'TCT': 0,
    'TGA': 0,
    'TGC': 0,
    'TGG': 0,
    'TGT': 0,
    'TTA': 0,
    'TTC': 0,
    'TTG': 0,
    'TTT': 0
}

# table for translating the codons into amino acids
table_of_codons = {
    'AAA': 'K',     # Lysine
    'AAC': 'N',     # Asparagine
    'AAG': 'K',     # Lysine
    'AAT': 'N',     # Asparagine
    'ACA': 'T',     # Threonine
    'ACC': 'T',     # Threonine
    'ACG': 'T',     # Threonine
    'ACT': 'T',     # Threonine
    'AGA': 'R',     # Arginine
    'AGC': 'S',     # Serine
    'AGG': 'R',     # Arginine
    'AGT': 'S',     # Serine
    'ATA': 'I',     # Isoleucine
    'ATC': 'I',     # Isoleucine
    'ATG': 'M',     # Methionine
    'ATT': 'I',     # Isoleucine
    'CAA': 'Q',     # Glutamine
    'CAC': 'H',     # Histidine
    'CAG': 'Q',     # Glutamine
    'CAT': 'H',     # Histidine
    'CCA': 'P',     # Proline
    'CCC': 'P',     # Proline
    'CCG': 'P',     # Proline
    'CCT': 'P',     # Proline
    'CGA': 'R',     # Arginine
    'CGC': 'R',     # Arginine
    'CGG': 'R',     # Arginine
    'CGT': 'R',     # Arginine
    'CTA': 'L',     # Leucine
    'CTC': 'L',     # Leucine
    'CTG': 'L',     # Leucine
    'CTT': 'L',     # Leucine
    'GAA': 'E',     # Glutamate
    'GAC': 'D',     # Aspartate
    'GAG': 'E',     # Glutamate
    'GAT': 'D',     # Aspartate
    'GCA': 'A',     # Alanine
    'GCC': 'A',     # Alanine
    'GCG': 'A',     # Alanine
    'GCT': 'A',     # Alanine
    'GGA': 'G',     # Glycine
    'GGC': 'G',     # Glycine
    'GGG': 'G',     # Glycine
    'GGT': 'G',     # Glycine
    'GTA': 'V',     # Valine
    'GTC': 'V',     # Valine
    'GTG': 'V',     # Valine
    'GTT': 'V',     # Valine
    'TAA': '*',     # Stop
    'TAC': 'Y',     # Tyrosine
    'TAG': '*',     # Stop
    'TAT': 'Y',     # Tyrosine
    'TCA': 'S',     # Serine
    'TCC': 'S',     # Serine
    'TCG': 'S',     # Serine
    'TCT': 'S',     # Serine
    'TGA': '*',     # Stop
    'TGC': 'C',     # Cystein
    'TGG': 'W',     # Tryptophane
    'TGT': 'C',     # Cysteine
    'TTA': 'L',     # Leucine
    'TTC': 'F',     # Phenylalanine
    'TTG': 'L',     # Leucine
    'TTT': 'F'     # Phenylalanine
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

# loop for checking every pair of start and stop
# "s" stands for single CDS
for s in segments_CDS:
    # "i" is an index according to whole sequence, in range from the start to the end of the CDS
    # range(start, stop, step)
    for i in range(s[0], s[1], 3):
        # building the codon = three subsequent nucleotides (string)
        single_codon = whole_sequence[i] + whole_sequence[i+1] + whole_sequence[i+2]
        # building the list of codons (string)
        codons_list.append(single_codon)

# print(codons_list)

# Counter builds a dictionary with number of occurrences as values
# .update() method allows to change values in the dictionary, leaving the unchanged ones
# if any triplet wasn't used in the given genome there will be corresponding value of 0
counted_codons.update(dict(Counter(codons_list)))

# print(counted_codons, '\n')

# now I need to reduce the number of Methionines by how many CDS I have = how many start codons
counted_codons['ATG'] -= counting_CDS
# print(counted_codons)

# adding the amino acid shortcut
# new dictionary with key = codon and value = list of two elements: amino acid shortcut and number of occurrences
for k in counted_codons.keys():
    for a in table_of_codons.keys():
        if k == a:
            temp_list = [table_of_codons[a], k, counted_codons[k]]
            list_of_counted_aminoacids.append(temp_list)
            # list_of_counted_aminoacids[k] = [table_of_codons[a], counted_codons[k]]

# print(list_of_counted_aminoacids)

with open('codon_usage.csv', 'w', newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter=',')
    my_writer.writerow(['aa', 'codon', 'occurrence'])
    my_writer.writerows(list_of_counted_aminoacids)

# testing my codons using original sequences from NCBI website
# two ranges for two CDS: 687-6057 or 5999-7928
test_list = []
# choose one of the below options to conduct a test for first or second CDS
try:
    test_no = int(input("choose test number: 1 or 2:\n"))
except ValueError:
    test_no = 0

# test no 1
if test_no == 1:
    for i in range(687, 6057, 3):
        # building the test codon
        test_codon = whole_sequence[i] + whole_sequence[i+1] + whole_sequence[i+2]
        for a in table_of_codons.keys():
            if test_codon == a:
                bla = [table_of_codons[a], test_codon, counted_codons[test_codon]]
                test_list.append(bla)
    # original sequence from NCBI website
    test_sequence = "MSGASSGTAIGAELFGISSVLGEYRVLIGDEGAGPSRSPSEVSFSVWYRSRAARLVVLCLVTSFLVPCLTFLIAEAVMGQTVTTPLSLTLDHWSEVRARAHNQGVEVRKKKWVTLCEAEWVMMNIGWPREGTFSLDNISQVEKKIFAPGPHGHPDQVPYITTWRSLATDPPSWVRPFLPPPKPPTPLPQPLSPQPSAPPTSSLYPVLPKPDPPKPPVLPPDPSSPLIDLLTEEPPPYPGGHGPPPSGPRTPAASPIASRLRERRENPAEESQALPLREGPNNRPQYWPFSASDLYNWKSHNPPFSQDPVALTNLIESILVTHQPTWDDCQQLLQALLTAEERQRVLLEARKQVPGEDGRPTQLPNVIDETFPLTRPNWDFATPAGREHLRLYRQLLLAGLRGAARRPTNLAQVKQVVQGKEETPAAFLERLKEAYRMYTPYDPEDPGQAASVILSFIYQSSPDIRNKLQRLEGLQGFTLSDLLKEAEKIYNKRETPEEREERLWQRQEERDKKRHKEMTKVLATVVAQNRDRDREESKLGDQRKIPLGKDQCAYCKEKGHWVRDCPKRPRKKPANSTLLNLEDXESQGQDPPPEPRITLRIGGQPVTFLVDTGAQHSVLTRPDGPLSDRTALVQGATGSRNYRWTTDRRVQLATGKETHSFLYVPECPYPLLGRDLLTKLKAQIHFTGEGANVVGPKGLPLQVLTLQLEEEYRLFEPESTLKQEMDIWLKNFPQAWAETGGIGMAHCQAPVLIQLKATATPISIRQYPMPHEAYQGIKPHIRRMLDQGILKPCRSPWNTPLLPVKKPGTEDYRPVQDLREVNKRVEDIHPTVPNPYNLLSTLPPSHPWYTVLDLKDAFFCLRLHPESQLLFAFEWKDPEIGLSGQLTWTRLPQGFKNSPTLFDEALHSDLADFRVRYPALVLLQYVDDLLLAAATKTECLEGTKALLETLGNKGYRASAKKAQICLQKVTYLGYSLEDGQRWLTKARKEAILSIPVPKNPRQVREFLGTAGYCRLWIPGFAELAAPLYPLTRPGTLFQWETKQQLAFENIKKALLSSPALGLPDITKPFELFIDENSGFAKGVLVQKLGPWKRPVAYLSKKLDTVASGWPPCLRMVAAIAILVKDADKLTLGQPLTILTSHPVEALVRQPPNKWLSNARMTHYQAMLLDAERVHFGPTVSLNPATLLPLPNEESHHDCLLVLAETHGTRPDLTDQPLPDADLTWYTDGSSFIRNGERKAGAAVTTESEVIWAASLPPGTSAQRAELIALTQALKMAKGKKLTVYTDSRYAFATAHVHGEIYRRRGLLTSEGKEIKNKNEILALLEALFLPKRLSIIHCPGHQKGDSPQAKGNRLADDTAKKAATETQSSLTILPTELIEGPKRPPWEYDNSDLDLVQKLEAHYEPKRGTWEYQGKTIMPEKYAKELISHLHKLTHLSARKMKTLLEREETGFYLPNRDLHLRQVTESCRACAQINAGKIKFGPDVRARGHRPGIHWEIDFTEIKPGMYGYKYLLVFIDTFSGWAEAYPAKHETAKVVAKKLLEEIFPRYGIPQVLGSDNGPAFISQVSQSVATLLGINWKLHCAYRPQSSGQVERMNRSIKETLTKLTLETGSKDWVLLLPLVLYRVRNTPGPHGLTPFEILYGAPPPLAHFFDADISSFATSPTMQAHLRALQLVQEEIQRPLAAAYREKLKTPVVPHPFKPGDSVWVRRHQTKNLEPRWKGPHIVLLTTPTALKVDGVAAWIHASHVKAAGPTTNQDPSDDPSSDDPSRWRVQRTQNPLKIRLSRGT"
    for t, s in zip(test_sequence, test_list):
        print(t, " -- ", s)
        print(t == s[0])
    # the NCBI sequence doesn't contain the stop triplet, so I have to print it additionally
    print(test_list[len(test_sequence)])

# test no 2
elif test_no == 2:
    for i in range(5999, 7928, 3):
        # building the test codon
        test_codon = whole_sequence[i] + whole_sequence[i+1] + whole_sequence[i+2]
        for a in table_of_codons.keys():
            if test_codon == a:
                bla = [table_of_codons[a], test_codon, counted_codons[test_codon]]
                test_list.append(bla)
    # original sequence from NCBI website
    test_sequence = "MEGSTHPKPSKDKTFSWDLMILVGVLLRLDVGMANPSPHQVHNVTWVITNVQKNSQANATSMLGTLTDAYPTLHVDLCDLVGDTWEPIVLNPNNVKHGARYSSSKYGCKTTDRKKQQQTYPFYVCPGHTPSMGPKGTHCGGAQDGFCAAWGCETTGEAWWKPTSSWDYITVKRGSSQDNSCEGKCNPLVLQFTQKGRQASWDGPKMWGLRLYRTGYDPVALFTVSRQVSTITPPQAMGPNLVLPDQKPPSRQSQTESKVATQKPQTNGSTPRSVAPATMGPKRIGTGDRLVNLVQGTYLALNATDPNKTRDCWLCLVSRPPYYEGIAILGNYSNQTNPPSSCLSTPQHKLTISEVSGQGLCIGTVPKTHQALCNKTQQGHTGAHYLAAPNGTYWACNTGLTPCISMAVLNWTSDFCVLIELWPRVTYHQPEYVYTHFDKAVRFRREPISLTVALMLGGLTVGGIAAGVGTGTKALLETAQFRQLQMAMHTDIQALEESISALERSLTSLSEVVLQNRRGLDILFLQEGGLCAALKEECCFYADHTGLVRDNMAKLRERLKQRQQLFDSQQGWFEGWFNKSPWFTTLISSIMGPLLILLLILLFGPCILNRLVQFVKDRISVVQALILTQQYQQIKQYDPDQP"
    for t, s in zip(test_sequence, test_list):
        print(t, " -- ", s)
        print(t == s[0])
    # the NCBI sequence doesn't contain the stop triplet, so I have to print it additionally
    print(test_list[len(test_sequence)])

else:
    pass

## AT4BR_Final_Project
### Kinga Zajdel

### Testing codon usage bias on the example of cat virus: FeLV (Feline Leukemia Virus)
Many of amino acids are coded by more than one nucleotide triplet. Silent mutations name such changes in DNA sequence that although altering a nucleotide(s) in the triplet, don't affect what amino acids are produced on the basis of the sequence. However, a change in the triplet may influence mRNA stability and thus alter number of protein molecules produced. Some findings suggest that those mutations aren't in fact silent and they can actually influence the gene expression.
I've chosen FeLV because it's a small genome, so data shouldn't load much computer's memory. My code should work also for bigger genomes.
In the middle of the first sequence there is a stop triplet (TAG). In some organisms it can encode Pyrrolysine but not in eukaryotes.

### How it works
Firstly, the Python script downloads the genome (sequence + additional data in Gene Bank format) from NCBI websitte, using given ID. According to provided CDS' start and end index nucleotide triplets are read and translated into amino acid abbreviations. Note that should any triplet be out of use in the given CDD it will get the value of 0. The data is exported into the .csv file that will be used by R Shiny Web App. Additionally, the script allows conducting a test which shows if amino acid sequence is accurate (according to data from NCBI).
Secondly, R Shiny Web App draws the data from the .csv file. User can choose amino acid which codon bias they want to see and kind of data representation: in number of occurrences in the genome or percentage of each triplet (100% is how many triplets were used to encode given amino acid). Adjusting the plot decoration is provided as well.

#### Sources
https://www.ncbi.nlm.nih.gov/nuccore/OR682571.1
"Coding-Sequence Determinants of Gene Expression in Escherichia coli", Grzegorz Kudla, Andrew W. Murray, David Tollervey, Joshua B. Plotkin
"The codon adaptation index - a measure of directional synonymous codon usage bias, and its potential applications", Paul M. Sharp, Wen-Hsiung Li
https://en.wikipedia.org/wiki/Pyrrolysine

/*No AI used./*

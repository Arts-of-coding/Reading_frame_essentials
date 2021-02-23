# Description Acseq (see markdown file)
The function Acseq generates a dictionary with keys as amino acids and values as codon triplets in three reading frames from input DNA sequences. The function generates a dictonary with amino acids and corresponding triplets. Additionally, the function generates a dataframe with sequences retaining the order of the original input, amino acids and frames. Please use DNA sequences with Uppercase "A", "C", "T" or "G". Other inputs will give a ValueError.

Functions determining the effect of point mutations and fusion proteins will be deposited in the future.

# Why choose Acseq?

There is no function that outputs the reading frames, sequences and codons together in one dataframe. Eventhough it is possible to get all these values from seperate occasions in Biopython, Acseq is more convenient. Although the Biopython's (https://biopython.org/docs/1.75/api/index.html) translation() function is versatile for different kinds of inputs, Acseq's translation part can also take sequences which are not multiples of 3. This entails that sequences do not need to be trimmed beforehand, pointing towards the simplicity and elegance of using Acseq. The aminoacid output of translate() from Biopyton has no joined values in the generated table (Acseqs adapted function Biopython function termed BIPY). Thus it is more difficult to use the output of the tabple in protein BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi). Furthermore, using the dictionary inside the function to generate the aminoacids is slightly faster for large datasets (>10000 entries with 1200 nucleotide length), see figure 1. In large appoaches, like RNA-seq or single cell RNA seq, datasets can be very large, so speed is an important aspect for analysing workflow.


# How to use the Acseq function?

First convert for example fastq files to fasta files with Biopython or another program in UNIX.
  for example:
  seqkit fq2fa reads_1.fq.gz -o reads_1.fa.gz #https://bioinf.shenwei.me/seqkit/usage/#fq2fa

Read in fasta files for generating the table: 
  from Bio import SeqIO
  seq_dict = {rec.id : rec.seq for rec in SeqIO.parse("myfile.fasta", "fasta")} # where myfile.fasta is your fasta file

From the dictionary with names and sequences you extract a list of only the sequences by:
  sequences = list(seq_dict.values())
    #note: this dictionary can also be used to find 
    # the SEQIDs in the fasta file corresponding to the sequences in the output table of Acseq

Now you can use the sequences as input for Acseq:
Acseq (sequences)

I recommend generating a dataframe out of the input function:
df = Acseq (sequences)

The dataframe can be if needed converted to an excel file:
df.to_excel("df.xlsx")  

# Help and support
The preferred way to get support is through the Github issues page.

Reach out to me at one of the following places!

GitHub
mail@julianarts.nl


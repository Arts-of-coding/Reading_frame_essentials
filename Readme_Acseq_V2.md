Functions determining the effect of point mutations and fusion proteins will be deposited in the future.

Description Acseq
The function Acseq generates a dictonary with amino acids and corresponding triplets. Please use DNA sequences with Uppercase "A", "C", "T" or "G". Other inputs will give a ValueError.

Why choose Acseq?
There is no function that outputs the reading frames, sequences and codons together in one dataframe. Eventhough it is possible to get all these values from seperate occasions in Biopython, Acseq is more convenient. Although the Biopython's (https://biopython.org/docs/1.75/api/index.html) translation() function is versatile for different kinds of inputs, Acseq's translation part can also take sequences which are not multiples of 3. This entails that sequences do not need to be trimmed beforehand, pointing towards the simplicity and elegance of using Acseq. The aminoacid output of translate() from Biopyton has no joined values in the generated table (Acseqs adapted function Biopython function termed BIPY). Thus it is more difficult to use the output of the tabple in protein BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi). Furthermore, using the dictionary inside the function to generate the aminoacids is slightly faster for large datasets, see benchmark figure. In large appoaches, like RNA-seq or single cell RNA seq, datasets can be very large, so speed is an important aspect for the general workflow.


How to use the Acseq function?
First convert for example fastq files to fasta files with Biopython or another program in UNIX.
  for example:
  seqkit fq2fa reads_1.fq.gz -o reads_1.fa.gz #https://bioinf.shenwei.me/seqkit/usage/#fq2fa

Read in fasta files for generating the table: 
  from Bio import SeqIO
  seq_dict = {rec.id : rec.seq for rec in SeqIO.parse("myfile.fasta", "fasta")} # where myfile.fasta is your fasta file

Or see example below
################################################################################################################################################################
Importing Fasta files the right way for analysis with Acseq


```python
from rfe import FastaRecord
from rfe import FastaParser
from rfe import lgene
from rfe import lID

sequence = lgene("YOURFILEHERE.fasta")
LIDs = lID("YOURFILEHERE.fasta")
```

Using Acseq to get the codons and corresponding aminoacids from the .fasta DNA sequences. First argument are the sequences and second argument are the IDs


```python
from rfe import Acseq
df = Acseq(sequence, LIDs)
```

Provided example


```python
from rfe import FastaRecord
from rfe import FastaParser
from rfe import lgene
from rfe import lID

sequence = lgene("test.fasta")
LIDs = lID("test.fasta")
```


```python
from rfe import Acseq
df = Acseq(sequence, LIDs)
```


```python
print (df[0:2])
```

                                             Amino_acids  Frame  \
    0  MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSS...      1   
    1  CLFFLFYCH*SLVSVLILQPELNYPLHTLILSHVVFITLTKFSDPQ...      2   
    
                                                  Codons  \
    0  [ATG, TTT, GTT, TTT, CTT, GTT, TTA, TTG, CCA, ...   
    1  [TGT, TTG, TTT, TTC, TTG, TTT, TAT, TGC, CAC, ...   
    
                                                ID#  
    0  >LC528232.1|surface glycoprotein|21566-25387  
    1  >LC528232.1|surface glycoprotein|21566-25387  
    

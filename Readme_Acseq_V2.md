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
    

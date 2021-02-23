# Acseq function
Version 4

Converts DNA sequences into amino acids and the corresponding codon triplets in three reading frames.


```python
def Acseq (seqs):

    # Importing nessesary libraries and specific functions
    import pandas as pd
    import re
    from re import match
    from pandas import DataFrame
    
    ### Set up for amino acids to codon convertion
    bases = ['T', 'C', 'A', 'G']
    codon = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codon, amino_acids))

    ### Creating empty list for frames function
    mycodons = []

    ### Making sure the input sequence is DNA
    for x in range(0,len(seqs)):
         for y in range(0,len(seqs[x])):
            if not match("^[ACTG]{1,}$", seqs[x][y],re.I):
                raise ValueError("List contains a RNA, non-DNA or typo in a sequence")
                
    ### Function runs if all sequences consist of DNA
    else:

    ### Determining for each of the three reading frames the triplets/codons
        def frames (seqs):
                codons = [seqs[i:i+3] for i in range(0,len(seqs),3)]
                mycodons.append(codons)
                codons2 = [seqs[j:j+3] for j in range(1,len(seqs),3)]
                mycodons.append(codons2)
                codons3 = [seqs[k:k+3] for k in range(2,len(seqs),3)]
                mycodons.append(codons3)
                    
        [frames (sequence) for sequence in seqs if sequence.upper().count('T') >= 0]
    
    ### Converting the triplets/codons to amino acids using the dictionary
        def replace_matched_items(word_list, dictionary):
            
            ### Omitting single and di-nucleotides from the conversion
            new_list = [[ele for ele in sub if ele != "T" if ele != "C" if ele != "A" 
                         if ele != "G" if ele != "TC" if ele != "CT" if ele != "CA" 
                         if ele != "AC" if ele != "TA" if ele != "AT" if ele != "CG" 
                         if ele != "GC" if ele != "TG" if ele != "GT" if ele != "TT" 
                         if ele != "TT" if ele != "CC" if ele != "AA" if ele != "GG"] for sub in word_list]
            ### 
            new_list = [[dictionary.get(item, item) for item in lst] for lst in new_list]
            return new_list
        
        triplet2 = replace_matched_items(mycodons, codon_table)
        myaa = [''.join(li) for li in triplet2]

        mySEQS = [[sequence for position in range(frame, len(sequence), 3)][0] for sequence in seqs for frame in range(3) if sequence.upper().count('T') >= 0]             
        myframes = [[frame + 1 for position in range(frame, len(sequence), 3)][0] for sequence in seqs for frame in range(3) if sequence.upper().count('T') >= 0]

    ### Generating readable names in the dataframe
        myframes2 = [" frame " + str(k) for k in myframes]
        res = [i + j for i, j in zip(mySEQS, myframes2)]
        
        Acseqdf = {'Sequence': mySEQS, 'Frame': myframes,'Amino_acids': myaa, 'Codons': mycodons}
        df = DataFrame(Acseqdf, columns = ['Amino_acids', 'Frame','Codons','Sequence'])
        return df

```


```python
import Bio.Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
```


```python
def BIPY (seqs):

    # Importing nessesary libraries and specific functions
    import pandas as pd
    import re
    from re import match
    from pandas import DataFrame
    
    ### Creating empty list for frames function
    mycodons = []

    ### Making sure the input sequence is DNA
    for x in range(0,len(seqs)):
         for y in range(0,len(seqs[x])):
            if not match("^[ACTG]{1,}$", seqs[x][y],re.I):
                raise ValueError("List contains a RNA, non-DNA or typo in a sequence")
                
    ### Function runs if all sequences consist of DNA
    else:

    ### Determining for each of the three reading frames the triplets/codons
        def frames (seqs):
                codons = [seqs[i:i+3] for i in range(0,len(seqs),3)]
                mycodons.append(codons)
                codons2 = [seqs[j:j+3] for j in range(1,len(seqs),3)]
                mycodons.append(codons2)
                codons3 = [seqs[k:k+3] for k in range(2,len(seqs),3)]
                mycodons.append(codons3)
                    
        [frames (sequence) for sequence in seqs if sequence.upper().count('T') >= 0]
    
    ### Bio python list comprehension for translation
    # this does not provide the right output, but computationally it is similar 
    # because the list has the same length
        myaa = [Seq(i).translate() for i in seqs for x in range(3)]
    
    ### List comprehensions for the original sequences and codons
        ###
        mySEQS = [[sequence for position in range(frame, len(sequence), 3)][0] for sequence in seqs for frame in range(3) if sequence.upper().count('T') >= 0]             
        myframes = [[frame + 1 for position in range(frame, len(sequence), 3)][0] for sequence in seqs for frame in range(3) if sequence.upper().count('T') >= 0]

    ### Generating readable names in the dataframe
        myframes2 = [" frame " + str(k) for k in myframes]
        res = [i + j for i, j in zip(mySEQS, myframes2)]
        
        Acseqdf = {'Sequence': mySEQS, 'Frame': myframes,'Amino_acids': myaa, 'Codons': mycodons}
        df = DataFrame(Acseqdf, columns = ['Amino_acids', 'Frame','Codons','Sequence'])
        return df
```


```python
# for res and res2
sequence = [None]*5000
import random

for n in range(5000):
    sequence[n] = ''
    for i in range(1500):
        sequence[n] += random.choice('ATCG')
```


```python
# https://stackoverflow.com/questions/33943362/timeit-equivalent-in-code
res = []
for i in range(3):
    a = %timeit -o Acseq(sequence)
    res.append(a)
```

    14.5 s ± 220 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    14.4 s ± 447 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    14.6 s ± 384 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    


```python
res2 = []
for i in range(3):
    a = %timeit -o BIPY(sequence)
    res2.append(a)
```

    17.8 s ± 820 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    17.3 s ± 592 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    17 s ± 514 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    


```python
res_best_times = [result.best for result in res] 
# time in seconds
print(res_best_times)
```

    [14.140888300000029, 13.718746099999976, 13.831045199999949]
    


```python
res2_best_times = [result.best for result in res2] 
# time in seconds
print(res2_best_times)
```

    [16.648856400000113, 16.300764899999876, 16.226890400000002]
    


```python
# for res3 and res4
sequence = [None]*15000
import random

for n in range(15000):
    sequence[n] = ''
    for i in range(1500):
        sequence[n] += random.choice('ATCG')
```


```python
res3 = []
for i in range(3):
    a = %timeit -o Acseq(sequence)
    res3.append(a)
```

    44 s ± 909 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    42.9 s ± 673 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    42.6 s ± 495 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    


```python
res4 = []
for i in range(3):
    a = %timeit -o BIPY(sequence)
    res4.append(a)
```

    48.9 s ± 538 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    47.7 s ± 278 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    50 s ± 1.07 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    


```python
res3_best_times = [result.best for result in res3] 
# time in seconds
print(res3_best_times)
```

    [42.63934029999973, 42.27674919999981, 41.972850200000266]
    


```python
res4_best_times = [result.best for result in res4] 
# time in seconds
print(res4_best_times)
```

    [48.23330370000076, 47.36441999999988, 48.34482440000011]
    


```python
# for res5 and res6
sequence = [None]*30000
import random

for n in range(30000):
    sequence[n] = ''
    for i in range(1500):
        sequence[n] += random.choice('ATCG')
```


```python
res5 = []
for i in range(3):
    a = %timeit -o Acseq(sequence)
    res5.append(a)
```

    1min 35s ± 2.67 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    1min 34s ± 1.71 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    1min 35s ± 1.79 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    


```python
res6 = []
for i in range(3):
    a = %timeit -o BIPY(sequence)
    res6.append(a)
```

    1min 47s ± 430 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    1min 47s ± 867 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    1min 44s ± 6.52 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    


```python
res5_best_times = [result.best for result in res5] 
# time in seconds
print(res5_best_times)
```

    [91.5960329999998, 92.3889063000006, 93.33658729999934]
    


```python
res6_best_times = [result.best for result in res6] 
# time in seconds
print(res6_best_times)
```

    [106.77714449999985, 106.23489049999989, 96.35497159999977]
    


```python
# Benchmarked on a i7 4790K processor
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
data1 = pd.DataFrame(np.random.rand(17,3), columns=['A','B','C']).assign(Location=1)
data2 = pd.DataFrame(np.random.rand(17,3)+0.2, columns=['A','B','C']).assign(Location=2)
data3 = pd.DataFrame(np.random.rand(17,3)+0.4, columns=['A','B','C']).assign(Location=3)
cdf = pd.concat([data1, data2, data3])
mdf = pd.melt(cdf, id_vars=['Location'], var_name=['Letter'])
ax = sns.barplot(x="Location", y="value", hue="Letter", data=mdf, errwidth=0)  
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, fancybox=True, shadow=True)
plt.show()
```

# Acseq function
Version 1.0

Converts DNA sequences into amino acids and the corresponding codon triplets in three reading frames.


```python
def Acseq (seqs):
    
    ### Set up for amino acids to codon convertion
    bases = ['T', 'C', 'A', 'G']
    codon = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codon, amino_acids))
    
    ### Creating empty lists for the dictionary
    myprotein = []
    mycodons = []
    
    ### Uncomment of code below is possible to see the codons and codon_table
    #print(codon)
    #print(codon_table)

    ### Making sure the input sequence is DNA
    count_a = 0
    count_c = 0
    count_g = 0
    count_t = 0
    
    for s in seqs:
        total = sum(len(s) for s in seqs)
        
    for x in range(0,len(seqs)):
         for y in range(0,len(seqs[x])):
            if "a" == seqs [x][y]:
                count_a = count_a + 1

    for a in range(0,len(seqs)):
        for b in range(0,len(seqs[a])):
            if "c"== seqs [a][b]:
                count_c = count_c + 1

    for c in range(0,len(seqs)):
        for d in range(0,len(seqs[c])):
            if "g" == seqs [c][d]:
                count_g = count_g + 1

    for e in range(0,len(seqs)):
        for f in range(0,len(seqs[e])):
            if "t" == seqs [e][f]:
                count_t = count_t + 1

    ### formula
    if total != count_a + count_c + count_g + count_t:
        
    ### can print the counts individually if you want to see them
        #print (count_a) 
        #print (count_c)
        #print (count_g)
        #print (count_t)
        #print (total)
        
        
        raise ValueError("List contains a RNA, non-DNA or typo in a sequence")
    
    ### Function runs if all sequences consist of DNA
    else:
        
    #Three nested functions for each of the codons in the reading frames 
        def frame_1 (seqs):
            codons = []
            for i in range(0,len(seqs),3): #3 for the length
            #print ("RF1 start base", i)       
                codons.append(seqs[i:i+3])
                #myDict["key2"].join(codons)
            print (codons)
            mycodons.append(codons)
            
        def frame_2 (seqs):
            codons2 = []
            for j in range(1,len(seqs),3): #3 for the length
            #print ("RF2 start base",j) 
                codons2.append(seqs[j:j+3])
                #myDict["key4"].append(codons2)
            print (codons2)
            mycodons.append(codons2)
            
        def frame_3 (seqs):
            codons3 = []
            for k in range(2,len(seqs),3): #3 for the length
            #print ("RF3 start base",k)    
                codons3.append(seqs[k:k+3])
                #myDict["key6"].append(codons3)
            print (codons3)
            mycodons.append(codons3)
            
    ### determining for each of the three reading frames the amino acids
        for sequence in seqs:
            if sequence.upper().count('T') >= 0 :
                for frame in range(3):
                    myDict = {} # Generating an empty dictionary
                    print ("\n"+"Reading frame " + str(frame+1 ) + " " + "&" + " " + "Start begin  of sequence: " + sequence [0:30]) 
                    protein = ''
                    for position in range(frame, len(sequence), 3):
                        triplet = sequence[position:position+3]
                        amino_acid = codon_table.get(triplet.upper(),'')
                        protein = protein+amino_acid
                    print(protein)
                    myprotein.append(protein)
                    #print(myprotein) see the list myprotein growing in each loop
                    
                    if (frame == 0):
                        frame_1(sequence)
                    elif (frame == 1):
                        frame_2(sequence)
                    elif (frame == 2):
                        frame_3(sequence)
    ### Filling the dictionary with the aminoacids as keys and the sequences as values
                    myDict = {key: value for key, value in zip(myprotein,mycodons)}
        return myDict

```


```python
### You can have your own list of sequences here

sequences = ['gatttcgggaattccggaattc','ccgccggaattcgaattc','attcgaccggaattcatgg'] #input list to see how the output is produced

Acseq(sequences)

```

    
    Reading frame 1 & Start begin  of sequence: gatttcgggaattccggaattc
    DFGNSGI
    ['gat', 'ttc', 'ggg', 'aat', 'tcc', 'gga', 'att', 'c']
    
    Reading frame 2 & Start begin  of sequence: gatttcgggaattccggaattc
    ISGIPEF
    ['att', 'tcg', 'gga', 'att', 'ccg', 'gaa', 'ttc']
    
    Reading frame 3 & Start begin  of sequence: gatttcgggaattccggaattc
    FREFRN
    ['ttt', 'cgg', 'gaa', 'ttc', 'cgg', 'aat', 'tc']
    
    Reading frame 1 & Start begin  of sequence: ccgccggaattcgaattc
    PPEFEF
    ['ccg', 'ccg', 'gaa', 'ttc', 'gaa', 'ttc']
    
    Reading frame 2 & Start begin  of sequence: ccgccggaattcgaattc
    RRNSN
    ['cgc', 'cgg', 'aat', 'tcg', 'aat', 'tc']
    
    Reading frame 3 & Start begin  of sequence: ccgccggaattcgaattc
    AGIRI
    ['gcc', 'gga', 'att', 'cga', 'att', 'c']
    
    Reading frame 1 & Start begin  of sequence: attcgaccggaattcatgg
    IRPEFM
    ['att', 'cga', 'ccg', 'gaa', 'ttc', 'atg', 'g']
    
    Reading frame 2 & Start begin  of sequence: attcgaccggaattcatgg
    FDRNSW
    ['ttc', 'gac', 'cgg', 'aat', 'tca', 'tgg']
    
    Reading frame 3 & Start begin  of sequence: attcgaccggaattcatgg
    STGIH
    ['tcg', 'acc', 'gga', 'att', 'cat', 'gg']
    




    {'DFGNSGI': ['gat', 'ttc', 'ggg', 'aat', 'tcc', 'gga', 'att', 'c'],
     'ISGIPEF': ['att', 'tcg', 'gga', 'att', 'ccg', 'gaa', 'ttc'],
     'FREFRN': ['ttt', 'cgg', 'gaa', 'ttc', 'cgg', 'aat', 'tc'],
     'PPEFEF': ['ccg', 'ccg', 'gaa', 'ttc', 'gaa', 'ttc'],
     'RRNSN': ['cgc', 'cgg', 'aat', 'tcg', 'aat', 'tc'],
     'AGIRI': ['gcc', 'gga', 'att', 'cga', 'att', 'c'],
     'IRPEFM': ['att', 'cga', 'ccg', 'gaa', 'ttc', 'atg', 'g'],
     'FDRNSW': ['ttc', 'gac', 'cgg', 'aat', 'tca', 'tgg'],
     'STGIH': ['tcg', 'acc', 'gga', 'att', 'cat', 'gg']}




```python
# You can make  a new dictionary with the function

new_dict = Acseq(sequences)
new_dict

```

    
    Reading frame 1 & Start begin  of sequence: gatttcgggaattccggaattc
    DFGNSGI
    ['gat', 'ttc', 'ggg', 'aat', 'tcc', 'gga', 'att', 'c']
    
    Reading frame 2 & Start begin  of sequence: gatttcgggaattccggaattc
    ISGIPEF
    ['att', 'tcg', 'gga', 'att', 'ccg', 'gaa', 'ttc']
    
    Reading frame 3 & Start begin  of sequence: gatttcgggaattccggaattc
    FREFRN
    ['ttt', 'cgg', 'gaa', 'ttc', 'cgg', 'aat', 'tc']
    
    Reading frame 1 & Start begin  of sequence: ccgccggaattcgaattc
    PPEFEF
    ['ccg', 'ccg', 'gaa', 'ttc', 'gaa', 'ttc']
    
    Reading frame 2 & Start begin  of sequence: ccgccggaattcgaattc
    RRNSN
    ['cgc', 'cgg', 'aat', 'tcg', 'aat', 'tc']
    
    Reading frame 3 & Start begin  of sequence: ccgccggaattcgaattc
    AGIRI
    ['gcc', 'gga', 'att', 'cga', 'att', 'c']
    
    Reading frame 1 & Start begin  of sequence: attcgaccggaattcatgg
    IRPEFM
    ['att', 'cga', 'ccg', 'gaa', 'ttc', 'atg', 'g']
    
    Reading frame 2 & Start begin  of sequence: attcgaccggaattcatgg
    FDRNSW
    ['ttc', 'gac', 'cgg', 'aat', 'tca', 'tgg']
    
    Reading frame 3 & Start begin  of sequence: attcgaccggaattcatgg
    STGIH
    ['tcg', 'acc', 'gga', 'att', 'cat', 'gg']
    




    {'DFGNSGI': ['gat', 'ttc', 'ggg', 'aat', 'tcc', 'gga', 'att', 'c'],
     'ISGIPEF': ['att', 'tcg', 'gga', 'att', 'ccg', 'gaa', 'ttc'],
     'FREFRN': ['ttt', 'cgg', 'gaa', 'ttc', 'cgg', 'aat', 'tc'],
     'PPEFEF': ['ccg', 'ccg', 'gaa', 'ttc', 'gaa', 'ttc'],
     'RRNSN': ['cgc', 'cgg', 'aat', 'tcg', 'aat', 'tc'],
     'AGIRI': ['gcc', 'gga', 'att', 'cga', 'att', 'c'],
     'IRPEFM': ['att', 'cga', 'ccg', 'gaa', 'ttc', 'atg', 'g'],
     'FDRNSW': ['ttc', 'gac', 'cgg', 'aat', 'tca', 'tgg'],
     'STGIH': ['tcg', 'acc', 'gga', 'att', 'cat', 'gg']}




```python

```

#!/usr/bin/env python
# coding: utf-8

# # Acseq function
# Version 1.0
# 
# Converts DNA sequences into amino acids and the corresponding codon triplets in three reading frames. Generates a dictionary and a dataframe.

# In[ ]:

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
            if not match("^[actg]{1,}$", seqs[x][y],re.I):
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
    
    ### Converting the codons to uppercase for conversion to amino acids using dictionary
        triplet2 = [[string.upper() for string in sublist] for sublist in mycodons]
        
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
        
        triplet2 = replace_matched_items(triplet2, codon_table)
        myaa = [''.join(li) for li in triplet2]
        
        mySEQS = [[sequence for position in range(frame, len(sequence), 3)][0] for sequence in seqs for frame in range(3) if sequence.upper().count('T') >= 0]             
        myframes = [[frame + 1 for position in range(frame, len(sequence), 3)][0] for sequence in seqs for frame in range(3) if sequence.upper().count('T') >= 0]

    ### Generating readable names in the dataframe
        Acseqdf = {'Sequence': mySEQS, 'Frame': myframes,'Amino_acids': myaa}
        df = DataFrame(Acseqdf, columns = ['Amino_acids', 'Frame','Sequence'])

    ### Generating the dictionary with sequences and corresponding codons          
        myDict = {key: value for key, value in zip(myaa,mycodons)}

        return myDict, df

# In[ ]:


### You can have your own list of sequences here

sequences = ['gatttcgggaattccggaattc','ccgccggaattcgaattc','attcgaccggaattcatgg'] #input list to see how the output is produced

Acseq(sequences)


# In[ ]:

# Generating a dataframe from the reading frames and sequences
new_df = Acseq(sequences)
new_df

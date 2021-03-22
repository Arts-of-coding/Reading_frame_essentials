class FastaRecord(object):
    """Class representing a FASTA record."""

    def __init__(self, description_line):
        """Initialise an instance of the FastaRecord class."""
        self.description = description_line.strip()
        self.sequences = []

    def add_sequence_line(self, sequence_line):
        """
        Add a sequence line to the FastaRecord instance.
        This function can be called more than once.
        """
        self.sequences.append( sequence_line.strip() )

    def matches(self, search_term):
        """Return True if the search_term is in the description."""
        return self.description.find(search_term) != -1

    def __repr__(self):
        """Representation of the FastaRecord instance."""
        lines = [self.description,]
        lines.extend(self.sequences)
        return '\n'.join(lines)

class FastaParser(object):
    """Class for parsing FASTA files."""

    def __init__(self, fpath):
        """Initialise an instance of the FastaParser."""
        self.fpath = fpath

    def __iter__(self):
        """Yield FastaRecord instances."""
        fasta_record = None
        with open(self.fpath, 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    if fasta_record:
                        yield fasta_record
                    fasta_record = FastaRecord(line)
                else:
                    fasta_record.add_sequence_line(line)
        yield fasta_record


def lgene(fyle):
    list_fasta = [''.join(fasta_record.sequences) for fasta_record in FastaParser(fyle) if not fasta_record.sequences == "Sequence unavailable"]
    return list_fasta

def lID(fyle):
    list_fasta2 = [fasta_record.description for fasta_record in FastaParser(fyle) if not fasta_record.sequences == "Sequence unavailable"]
    return list_fasta2

def Acseq(seqs, LIDs):

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
    NLIDs = [item for item in LIDs for i in range(3)]

    ### Making sure the input sequence is DNA
    for x in range(0,len(seqs)):
         for y in range(0,len(seqs[x])):
            if not match("^[ACTGNBDHKMNRSVWY]{1,}$", seqs[x][y],re.I):
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
        
        myframes = [[frame + 1 for position in range(frame, len(sequence), 3)][0] for sequence in seqs for frame in range(3) if sequence.upper().count('T') >= 0]
    ### Checking if fasta.py is used for changing of the identifiers
            
    ### Generating readable names in the dataframe
        myframes2 = [" frame " + str(k) for k in myframes]
        res = [i + j for i, j in zip(LIDs, myframes2)]
        
        Acseqdf = {'ID#': NLIDs, 'Frame': myframes,'Amino_acids': myaa, 'Codons': mycodons}
        df = DataFrame(Acseqdf, columns = ['Amino_acids', 'Frame','Codons','ID#'])
        return df
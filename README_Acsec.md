

### or if you want to do it manually with: #source: https://www.tjelvarolsson.com/blog/object-oriented-programming-for-scientists/

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

    def __repr__(self):
        """Representation of the FastaRecord instance."""
        lines = [self.description,]
        lines.extend(self.sequences)
        return '\n'.join(lines)
        
    
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


# Syntax requires python version >= 3.6.x

from pathlib import Path


class FASTA:
    """FASTA file parser
    Attrs:
        file_path (str): path to the FASTA file

        sequence (Sequence): nucleotide sequence

        version (str): nucleotide sequence identification 
        number that represents a single, specific 
        sequence in the GenBank database

        definition (str): brief description of sequence, includes 
        information such as source organism, gene name/protein name, 
        or some description of the sequence's function
    """

    def __init__(self, file_path: str):
        """
        Args:
            file_path (str): Path to FASFA file
        """

        self.file_path = Path(file_path).resolve()
        self.sequence: Sequence = None
        self.version: str = ""
        self.definition: str = ""

        self.parse_fasta_file()

    def __repr__(self):
        """String representation of FASTA object

        Used automatically when printing the object

        Returns:
            str : string representation of FASTA object
        """

        return "<{}> {}".format(self.version, self.definition)

    def parse_fasta_file(self):
        """Parse FASTA file

        Should not be invoked directly

        Note: The sequence does not have to be reversed because
        all FASTA files read from 5' to 3'
        """

        file = open(self.file_path, "r")
        header = file.readline().rstrip('\n')

        end_of_version = header.index(" ")
        self.version = header[1:end_of_version]
        self.definition = header[end_of_version+1:]

        bases = ""
        line = file.readline().rstrip('\n')
        while(line):
            bases += line
            line = file.readline().rstrip('\n')

        self.sequence = Sequence(bases)

        file.close()


class Sequence:
    """A sequence of nucleotides"""

    def __init__(self, bases: str):
        self.bases: str = bases or ""

    def __repr__(self):
        """String representation of Sequence object

        Used automatically when printing the object
        Truncates the bases to the first 10 with an ellipses
        if the length exceeds 10

        Returns:
            str : string representation of Sequence object
        """

        preview_len = 10
        return "<Sequence> 5'-{}{}-3'".format(
            self.bases[0:preview_len],
            "..." if len(self.bases) > preview_len else "")

    def slice(self, begin: int, end: int = None):
        """Create a subsequence of nucleotides

        Begin and end indices start at 1 to be consistent with GenBank

        Args:
            begin (int): the position of the start of the subsequence
            end (int, optional): the position of the end of the subsequence. Defaults to None.

        Returns:
            Sequence: A subsequence of the original Sequence object
        """
        return Sequence(self.bases[begin-1:end])

    def compliment(self):
        """Create a compliment of the sequence

        Compliment has the complimentary base pairs and is flipped to maintain
        a reading of 5' to 3'
        """
        c_strand = ""

        for base in self.bases:
            if base == 'A':
                c_strand += 'T'
            if base == 'T':
                c_strand += 'A'
            if base == 'C':
                c_strand += 'G'
            if base == 'G':
                c_strand += 'C'

        return Sequence(c_strand[::-1])

    def gc_content(self):
        """Computes the GC content of a sequence

        GC content is the ratio of G and C bases to the total sequence length
        """
        return (self.bases.count('C') + self.bases.count('G')) / len(self.bases)

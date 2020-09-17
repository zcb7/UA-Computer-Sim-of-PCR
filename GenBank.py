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
        return "<Sequence> {}{}".format(
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

    def group_codons(self):
        nucleotides_in_codon = 3
        codon_sequences = [self.bases[i:i+nucleotides_in_codon]
                           for i in range(0, len(self.bases), nucleotides_in_codon)]
        return [Codon.from_sequence(x) for x in codon_sequences]

    def compliment(self):
        c_strand = ""

        for base in self.bases:
            if base is 'A':
                c_strand += 'T'
            if base is 'T':
                c_strand += 'A'
            if base is 'C':
                c_strand += 'G'
            if base is 'G':
                c_strand += 'C'

        return Sequence(c_strand)


class Codon:
    def __init__(self, name: str, abbreviation: str, symbol: str, start=False, stop=False):
        self.name = name
        self.abbreviation = abbreviation
        self.symbol = symbol
        self.start = start
        self.stop = stop

    def __repr__(self):
        return "({}/{}) {}".format(self.abbreviation, self.symbol, self.name)

    @staticmethod
    def from_sequence(sequence: str):
        if(len(sequence) is not 3):
            raise 'Sequence must be of length 3'

        return CODON_LOOKUP_TREE[sequence[0]][sequence[1]][sequence[2]]


CODON_LOOKUP_TREE = {
    "A": {
        "A": {
            "A": Codon('Lysine', 'Lys', 'K'),
            "T": Codon('Asparagine', 'Asn', 'N'),
            "C": Codon('Asparagine', 'Asn', 'N'),
            "G": Codon('Lysine', 'Lys', 'K')
        },
        "T": {
            "A": Codon('Isoleucine', 'Ile', 'I'),
            "T": Codon('Isoleucine', 'Ile', 'I'),
            "C": Codon('Isoleucine', 'Ile', 'I'),
            "G": Codon('Methionine', 'Met', 'M', start=True)
        },
        "C": {
            "A": Codon('Threonine', 'Thr', 'T'),
            "T": Codon('Threonine', 'Thr', 'T'),
            "C": Codon('Threonine', 'Thr', 'T'),
            "G": Codon('Threonine', 'Thr', 'T')
        },
        "G": {
            "A": Codon('Arginine', 'Arg', 'R'),
            "T": Codon('Serine', 'Ser', 'S'),
            "C": Codon('Serine', 'Ser', 'S'),
            "G": Codon('Arginine', 'Arg', 'R')
        }
    },
    "T": {
        "A": {
            "A": Codon('Ochre', None, None, stop=True),
            "T": Codon('Tyrosine', 'Tyr', 'Y'),
            "C": Codon('Tyrosine', 'Tyr', 'Y'),
            "G": Codon('Amber', None, None, stop=True)
        },
        "T": {
            "A": Codon('Leucine', 'Leu', 'L'),
            "T": Codon('Phenylalanine', 'Phe', 'F'),
            "C": Codon('Phenylalanine', 'Phe', 'F'),
            "G": Codon('Leucine', 'Leu', 'L')
        },
        "C": {
            "A": Codon('Serine', 'Ser', 'S'),
            "T": Codon('Serine', 'Ser', 'S'),
            "C": Codon('Serine', 'Ser', 'S'),
            "G": Codon('Serine', 'Ser', 'S')
        },
        "G": {
            "A": Codon('Opal', None, None, stop=True),
            "T": Codon('Cysteine', 'Cys', 'C'),
            "C": Codon('Cysteine', 'Cys', 'C'),
            "G": Codon('Tryptophan', 'Trp', 'W')
        }
    },
    "C": {
        "A": {
            "A": Codon('Glutamine', 'Gln', 'Q'),
            "T": Codon('Histidine', 'His', 'H'),
            "C": Codon('Histidine', 'His', 'H'),
            "G": Codon('Glutamine', 'Gln', 'Q')
        },
        "T": {
            "A": Codon('Leucine', 'Leu', 'L'),
            "T": Codon('Leucine', 'Leu', 'L'),
            "C": Codon('Leucine', 'Leu', 'L'),
            "G": Codon('Leucine', 'Leu', 'L')
        },
        "C": {
            "A": Codon('Proline', 'Pro', 'P'),
            "T": Codon('Proline', 'Pro', 'P'),
            "C": Codon('Proline', 'Pro', 'P'),
            "G": Codon('Proline', 'Pro', 'P')
        },
        "G": {
            "A": Codon('Arginine', 'Arg', 'R'),
            "T": Codon('Arginine', 'Arg', 'R'),
            "C": Codon('Arginine', 'Arg', 'R'),
            "G": Codon('Arginine', 'Arg', 'R')
        }
    },
    "G": {
        "A": {
            "A": Codon('Glutamic acid', 'Glu', 'E'),
            "T": Codon('Aspartic acid', 'Asp', 'D'),
            "C": Codon('Aspartic acid', 'Asp', 'D'),
            "G": Codon('Glutamic acid', 'Glu', 'E')
        },
        "T": {
            "A": Codon('Valine', 'Val', 'V'),
            "T": Codon('Valine', 'Val', 'V'),
            "C": Codon('Valine', 'Val', 'V'),
            "G": Codon('Valine', 'Val', 'V')
        },
        "C": {
            "A": Codon('Alanine', 'Ala', 'A'),
            "T": Codon('Alanine', 'Ala', 'A'),
            "C": Codon('Alanine', 'Ala', 'A'),
            "G": Codon('Alanine', 'Ala', 'A')
        },
        "G": {
            "A": Codon('Glycine', 'Gly', 'G'),
            "T": Codon('Glycine', 'Gly', 'G'),
            "C": Codon('Glycine', 'Gly', 'G'),
            "G": Codon('Glycine', 'Gly', 'G')
        }
    }
}

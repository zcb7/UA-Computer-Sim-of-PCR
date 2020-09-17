from GenBank import FASTA


def main():
    genome = FASTA('./data/Wuhan-Hu-1_complete_genome.fasta')

    # Indices from GenBank for M gene
    M_gene = genome.sequence.slice(26523, 27191)

    print('M Gene Sequence\n{}\n'.format(M_gene.bases))

    # I think this needs to be flipped to keep all sequences in the same direction
    # (5' -> 3' / 3' -> 5')
    print('M Gene Sequence cDNA\n{}\n'.format(M_gene.compliment().bases))

    # I don't think we will be using this part in the project
    symbols = map(lambda codon: codon.symbol or "", M_gene.group_codons())
    protein_string = "".join(symbols)
    print('Protein String\n{}\n'.format(protein_string))


if __name__ == "__main__":
    main()

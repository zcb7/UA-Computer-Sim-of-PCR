from GenBank import FASTA


def main():
    genome = FASTA('./data/Wuhan-Hu-1_complete_genome.fasta')

    # Indices from GenBank for M gene
    M_gene = genome.sequence.slice(26523, 27191)

    print('M Gene Sequence\n{}\n'.format(M_gene.bases))

    print('M Gene Sequence cDNA\n{}\n'.format(M_gene.compliment().bases))

    print('GC%: {}'.format(M_gene.gc_content()))


if __name__ == "__main__":
    main()

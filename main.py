from time import time

from GenBank import FASTA, Sequence
from PCR import PCR


def main():
    genome = FASTA('./data/Wuhan-Hu-1_complete_genome.fasta')

    # Indices from GenBank for M gene
    M_gene = genome.sequence.slice(26523, 27191)

    # Creates tuple RNA/cDNA tuple
    paired_sequence = (M_gene, M_gene.compliment())

    # Primer sequences, primer pair 3 from:
    # https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1600372805&job_key=iIJXWm09YJVHq3qud85enA3VT64gxlSzIQ&CheckStatus=Check
    primers = (Sequence('GCTTGTTTTGTGCTTGCTGC'),
               Sequence('GGAGTGGCACGTTGAGAAGA'))

    # Other PCR parameters
    num_cycles = 20
    fall_off_noise = 50

    # Begin PCR
    start_time = time()
    results = PCR(paired_sequence, primers, num_cycles=num_cycles,
                  fall_off_noise=fall_off_noise)
    end_time = time()
    time_elapsed = round(end_time - start_time, 2)
    print('Completed {} total cycles in {} seconds'.format(
        num_cycles, time_elapsed))

    # Write results to file
    out_file = open('result.csv', 'w')
    out_file.write('strand,complimentary strand\n')
    for result in results:
        if result[0] != None:
            out_file.write("{}".format(result[0].bases))
        if result[1] != None:
            out_file.write("{}".format(result[1].bases))
        out_file.write('\n')
    out_file.close()
    print('Results written to result.csv')

    # TODO Plot results with matplotlib and interpret


if __name__ == "__main__":
    main()

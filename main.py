from GenBank import FASTA

# def PCR(paired_sequence, f_primer, r_primer, cycles):

def main():
    genome = FASTA('./data/Wuhan-Hu-1_complete_genome.fasta')

    # Indices from GenBank for M gene
    M_gene = genome.sequence.slice(26523, 27191)

    #Creates tuple RNA/cDNA tuple
    paired_sequence = [[M_gene, M_gene.compliment()]]

    #Primer pair 3 from
    #https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1600372805&job_key=iIJXWm09YJVHq3qud85enA3VT64gxlSzIQ&CheckStatus=Check
    
    #Primer sequence, Start point, Stop point
    f_primer = ("GCTTGTTTTGTGCTTGCTGC") #, 187, 206)
    f_primer_comp = "CGAACAAAACACGAACGACG"

    r_primer = ("GGAGTGGCACGTTGAGAAGA") #, 373, 354)
    r_primer_comp = "CCTCACCGTGCAACTCTTCT"
    r_primer_rev_comp = "TCTTCTCAACGTGCCACTCC"

    cycles = 20
    i = 0

    while i < 2 ** cycles:
        if i == 0:
            if paired_sequence[i][0].bases.find(f_primer) == -1 or paired_sequence[i][1].bases.find(r_primer) == -1:
                break
            else:
                f_start = paired_sequence[i][0].bases.find(f_primer) + 1
                r_start = paired_sequence[i][1].bases.find(r_primer) + 1
                
                print(f_start)
                print(r_start)

                f_comp = paired_sequence[i][0].slice(f_start,).compliment()
                r_comp = paired_sequence[i][1].slice(r_start,).compliment()

                print(paired_sequence[i][0].slice(f_start,).bases)
                print('\n')
                print(f_comp.bases)
                print('\n')
                print(paired_sequence[i][1].slice(r_start,).bases)
                print('\n')
                print(r_comp.bases)
                print('\n')

                paired_sequence.append([paired_sequence[i][0].slice(f_start,), paired_sequence[i][0].slice(f_start,).compliment()])
                paired_sequence.append([paired_sequence[i][1].slice(r_start,), paired_sequence[i][1].slice(r_start,).compliment()])
                
                i += 1
                
        elif paired_sequence[i][0].bases.find(f_primer) != -1:

            f_start = paired_sequence[i][0].bases.find(f_primer) + 1
            r_start = paired_sequence[i][1].bases.find(r_primer) + 1
            
            print(r_start)

            r_comp = paired_sequence[i][1].slice(r_start,).compliment()

            print(paired_sequence[i][1].slice(r_start,).bases)
            print('\n')
            print(r_comp.bases)
            print('\n')

            paired_sequence.append([paired_sequence[i][1].slice(r_start,), paired_sequence[i][1].slice(r_start,).compliment()])

            i += 1

        else:
            if paired_sequence[i][1].bases.find(f_primer) == -1 or paired_sequence[i][0].bases.find(r_primer) == -1:
                break
            else:
                f_start = paired_sequence[i][1].bases.find(f_primer) + 1
                r_start = paired_sequence[i][1].bases.find(r_primer) + 1
                
                print(f_start)

                f_comp = paired_sequence[i][1].slice(f_start,).compliment()

                print(paired_sequence[i][1].slice(f_start,).bases)
                print('\n')
                print(f_comp.bases)
                print('\n')

                paired_sequence.append([paired_sequence[i][1].slice(f_start,), paired_sequence[i][1].slice(f_start,).compliment()])

                i += 1

    print('M Gene Sequence\n{}\n'.format(M_gene.bases))

    print('M Gene Sequence cDNA\n{}\n'.format(M_gene.compliment().bases))

    print('GC%: {}'.format(M_gene.gc_content()))


if __name__ == "__main__":
    main()

    

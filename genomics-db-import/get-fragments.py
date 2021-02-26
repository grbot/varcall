#!/usr/bin/env python3

from optparse import OptionParser
import csv, sys, os, math

# genome-size-file can be created like this
# cat /cbio/dbs/gatk/2.8/b37/human_g1k_v37_decoy.fasta.fai | cut -f 1,2 | grep -P "^1\t|^2\t|^3\t|^4\t|^5\t|^6\t|^7\t|^8\t|^9\t|^10\t|^11\t|^12\t|^13\t|^14\t|^15\t|^16\t|^17\t|^18\t|^19\t|^20\t|^21\t|^22\t|^X\t|^Y\t|^MT\t" > genome.sizes

def main():
    parser = OptionParser(usage="usage: %prog -s sample-size -g genome-size-file -i intervals-file")
    parser.add_option("-s", "--sample-size", dest="sample_size", help="Total size of samples in GenomicsDBImport.")
    parser.add_option("-g", "--genome-size-file", dest="genome_size_file", help="Chromosomes and corresponding sizes (e.g. 1\t12312123).")
    parser.add_option("-i", "--intervals-file", dest="intervals_file", help="Output intervals files to be used as input to GenomicsDBImport.")
    (options, args) = parser.parse_args()
    if not options.sample_size:
        print("No sample size specified, please specify with the -s flag.")
        return -1
    if not options.genome_size_file:
        print("No genome size file specified, please specify with the -g flag.")
        return -2
    if not options.intervals_file:
        print("No intervals file specified, please specify with the -i flag.")
        return -3
    
    sample_size = (int)(options.sample_size)
    genome_size_file = options.genome_size_file
    intervals_file = options.intervals_file

    genome_sizes = csv.reader(open(genome_size_file), delimiter='\t')
    intervals = open(intervals_file, 'w')

    total = 0
    genome = {}
    for row in genome_sizes:
        total+=int(row[1])
        genome[row[0]] = row[1]

    fragment_size = math.floor(total / sample_size)

    for chr in genome:
        interval_start = 1
        interval_end = interval_start + fragment_size
        while interval_end <= int(genome[chr]):
            intervals.write (chr + ":" + str(interval_start) + "-" + str(interval_end) + "\n")
            interval_start = interval_end + 1
            interval_end = interval_start + fragment_size
        if (interval_end != int(genome[chr])):
           interval_end = int(genome[chr])
           intervals.write (chr + ":" + str(interval_start) + "-" + str(interval_end) + "\n")

    intervals.close()

if __name__ == "__main__":
    sys.exit(main())

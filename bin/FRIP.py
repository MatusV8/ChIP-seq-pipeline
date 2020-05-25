#!/usr/bin/env python

# import modules
import deeptools.countReadsPerBin as crpb
import pysam
import argparse

# Load input files
parser = argparse.ArgumentParser(description = 'This is a script to calculate % of reads in a bed file')
parser.add_argument("--bam", help = "input bam file")
parser.add_argument("--bed", help = "input bed file")
parser.add_argument("-n", "--name", help = "sample name")
parser.add_argument("-p", "--processors", help = "number of processors", type = int)
args = parser.parse_args()

# Do not calculate if the bedfile is empty

num_lines = 0
with open(args.bed, 'r') as f:
    for line in f:
        num_lines += 1

# Calculate Reads in bam file
if num_lines > 0:
    cr = crpb.CountReadsPerBin([args.bam], bedFile= args.bed, numberOfProcessors= args.processors)
    reads_at_peaks = cr.run()

    # Calculate total number of reads in peaks
    total = reads_at_peaks.sum(axis=0)

    # Calculate % of frangments in peaks
    bam = pysam.AlignmentFile(args.bam)
    frip = float(total[0]) / bam.mapped


    print(str(frip*100))
else:
    print('0')

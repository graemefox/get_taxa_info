#!/usr/bin/python -tt

#### remove lines from a fasta file that do not have any sequence (empty ID lines)

import argparse
parser = argparse.ArgumentParser(description='pal_filter')
parser.add_argument('-i','--input1', help='taxonomy file', \
                                    required=True)
parser.add_argument('-o','--output', help='biom file in txt format', \
                                    required=True)
args = parser.parse_args()
IDs = []
Sequences = []
with open(args.input1) as input_file:
    for line in input_file:
        if line.startswith(">"):
            IDs.append(line)
        else:
            Sequences.append(line)

with open(args.output, 'w') as output_file:
    for x, y in zip(IDs, Sequences):
        if len(y) > 1:
            output_file.write(x + y)

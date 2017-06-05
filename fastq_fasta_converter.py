#!/usr/bin/python -tt

import argparse, os, Bio
from Bio import SeqIO

parser = argparse.ArgumentParser(description='pal_filter')
parser.add_argument('-i','--input1', help='Fastq File', \
                                    required=True)
parser.add_argument('-o','--input2', help='output file', \
                                    required=True)

args = parser.parse_args()

sample_name = (args.input1[0:18])
#print(sample_name)


#with open (args.input1) as file:
with open(args.input1) as f:
    content = f.readlines()

no_plus = []

for entry in content:
    if len(entry.rstrip("\n")) != 1:
        no_plus.append(entry)

with open (args.input2, 'w') as out:
    count = 0
    while count < len(no_plus):
        #out.write(no_plus[count].replace("@HISEQ", ">HISEQ"))
        out.write(">" + sample_name + str(count) +"\n")
#        print(no_plus[count].replace("@HISEQ", ">HISEQ"))
        out.write(no_plus[count+1])
        count = count + 3

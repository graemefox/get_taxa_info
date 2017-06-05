#!/usr/bin/python -tt

import re, argparse

# something to try to sort out the reference database and taxonomy files

parser = argparse.ArgumentParser(description='pal_filter')
parser.add_argument('-i','--input1', help='Fastq File', \
                                    required=True)
parser.add_argument('-o','--input2', help='sequence output file', \
                                    required=True)
parser.add_argument('-p','--input3', help='taxonomy output file', \
                                    required=True)

args = parser.parse_args()

with open(args.input1, 'r') as input_file:
    content = input_file.readlines()

with open(args.input2, 'w') as outputfile:
    with open(args.input3, 'w') as tax_out_file:
        for row in content:

            if row.startswith(">"):
                count = 0
                matchObj = re.match( r'>[0-9]*;', row)
                if matchObj:
                    outputfile.write(matchObj.group(0).rstrip(";") + "\n")
                while count < len(row.split(";")):
                    if count == 0:
                        tax_out_file.write(row.split(";")[count].lstrip(">") + " ")
                    if not(count == 1 or count == 2):
                        data = row.split(";")[count].lstrip(">")
                        matchObj = re.match( r'[a-z]__[A-Za-z]*', data)
                        if matchObj:
                            tax_out_file.write((matchObj.group(0))+"; ")
                    count = count + 1
                tax_out_file.write("\n")
            else:
                outputfile.write(row)

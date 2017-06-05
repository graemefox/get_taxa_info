#!/usr/bin/python -tt
#FUNCTIONS

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

import argparse
import Bio
from collections import defaultdict
import os.path
from Bio import SeqIO
import subprocess
import re
# parse input arguments
parser = argparse.ArgumentParser(description='pal_filter')
parser.add_argument('-i','--input1', help='taxonomy file', \
                                    required=True)
parser.add_argument('-j','--input2', help='biom file in txt format', \
                                    required=True)
parser.add_argument('-k','--input3', help='new_ref_seqs.fna file associated with this Sample', \
                                    required=True)
parser.add_argument('-r','--input4', help='sequence file of the taxonomy database', \
                                    required=True)
parser.add_argument('-o','--output', help='output file containing taxa info', \
                                    required=True)
parser.add_argument('-p','--summary_output', help='file containing summary of taxa info', \
                                    required=True)
parser.add_argument('-s','--size_of_ref_fragment', help='size of total barcoding region (not your amplicon)/ IE. 16S is 1500.', \
                                    required=True)
parser.add_argument('-t','--taxidnucl', help='file of taxa info provided by NCBI (ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)', \
                                    required=True)
parser.add_argument('-u','--names', help='names.dmp provided by NCBI (ftp://ftp.ncbi.nih.gov/pub/taxonomy/)', \
                                    required=True)
parser.add_argument('-v','--nodes', help='nodes.dmp provided by NCBI (ftp://ftp.ncbi.nih.gov/pub/taxonomy/)', \
                                    required=True)
args = parser.parse_args()

# read in everything we need from the input files so that it can be easily
# written to the output
# get all the taxa info into a dicts
dicts = []
with open(args.input1, 'r') as input_file:
    for line in input_file:
        dicts.append({'ID': line.split("\t")[0], 'Taxa': line.split("\t")[1]})
input_file.close()

# get all new_ref_seqs.fna into a dictionary
ref_seqs_dict = SeqIO.to_dict(SeqIO.parse(args.input3, "fasta"))
taxonomy_seq = SeqIO.to_dict(SeqIO.parse(args.input4, "fasta"))

# get total seq count:
total = 0
with open(args.input2, 'r') as count_total:
    next(count_total)
    next(count_total)
    for line in count_total:
        total = total + float(line.split("\t")[1])
total = (int(total))

# define some variables for the blast searches
program = "blastn"
database = "nt"
E_VALUE_THRESH = 0.04

# main grabbing of info and writing of summary output file:
KINGDOM_COUNTER = {}
PHYLA_COUNTER = {}
CLASS_COUNTER = {}
ORDER_COUNTER = {}
FAMILY_COUNTER = {}
GENUS_COUNTER = {}
SPECIES_COUNTER = {}

with open(args.output, 'w') as output_file:
    output_file.write("OTU_ID\tSequence_Counts\tRelative_Abundance(%)\tTaxonomy\tReference_Sequence\n")
    list_of_data_to_blast = []
    total_lengths = 0
    list_of_fragment_sizes = []
    with open("temp_blast_input.fasta", 'w') as blast_input:
        with open(args.input2, 'r') as biom_file:
            next(biom_file)
            next(biom_file)
            for line in biom_file:
                print(line)
                success = "0"
                for row in dicts:
                    if row['ID'] == line.split("\t")[0]:
                        # for each line in the file it now attempts to get taxonomy info in 3 different ways:
                        ##################################
                        # OPTION 1  -  confidence level #1
                        ##################################
                        # if it has a taxonomy entry in the provided taxonomy file, then grab that:
                        success = "1"
                        rel_abund = 100*(float(line.split("\t")[1])/total)
                        seq = taxonomy_seq[line.split("\t")[0]]
                        sequence = seq.seq
                        # remove any stray newlines from the taxa line:
                        row['Taxa'] = row['Taxa'].replace("\n", "")
                        output_file.write(line.split("\t")[0] + "\t" + line.split("\t")[1].rstrip("\n") + "\t" + str(rel_abund) + "\t" + row['Taxa'].rstrip("\n") + "\t" + str(sequence) + "\n")
                        # the taxa info from the provided taxa file is then used to construct
                        #the summary stats as this is the stuff we are most confident in
                        # do all the taxa counting. tedious looping but it runs.
                        # if it does already have an entry in the counter, create one
                        if not row['Taxa'].split("; ")[0].lstrip("k__") in KINGDOM_COUNTER:
                            KINGDOM_COUNTER[row['Taxa'].split("; ")[0].lstrip("k__")] = 1
                        else:
                            KINGDOM_COUNTER[row['Taxa'].split("; ")[0].lstrip("k__")] += 1
                        if not row['Taxa'].split("; ")[1].lstrip("p__") in PHYLA_COUNTER:
                            PHYLA_COUNTER[row['Taxa'].split("; ")[1].lstrip("p__")] = 1
                        else:
                            PHYLA_COUNTER[row['Taxa'].split("; ")[1].lstrip("p__")] += 1
                        if not row['Taxa'].split("; ")[2].lstrip("c__") in CLASS_COUNTER:
                            CLASS_COUNTER[row['Taxa'].split("; ")[2].lstrip("c__")] = 1
                        else:
                            CLASS_COUNTER[row['Taxa'].split("; ")[2].lstrip("c__")] += 1
                        if not row['Taxa'].split("; ")[3].lstrip("o__") in ORDER_COUNTER:
                            ORDER_COUNTER[row['Taxa'].split("; ")[3].lstrip("o__")] = 1
                        else:
                            ORDER_COUNTER[row['Taxa'].split("; ")[3].lstrip("o__")] += 1
                        if not row['Taxa'].split("; ")[4].lstrip("f__") in FAMILY_COUNTER:
                            FAMILY_COUNTER[row['Taxa'].split("; ")[4].lstrip("f__")] = 1
                        else:
                            FAMILY_COUNTER[row['Taxa'].split("; ")[4].lstrip("f__")] += 1
                        if not row['Taxa'].split("; ")[5].lstrip("g__") in GENUS_COUNTER:
                            GENUS_COUNTER[row['Taxa'].split("; ")[5].lstrip("g__")] = 1
                        else:
                            GENUS_COUNTER[row['Taxa'].split("; ")[5].lstrip("g__")] += 1
                        if not row['Taxa'].split("; ")[6].lstrip("s__") in SPECIES_COUNTER:
                            SPECIES_COUNTER[row['Taxa'].split("; ")[6].lstrip("s__")] = 1
                        else:
                            SPECIES_COUNTER[row['Taxa'].split("; ")[6].lstrip("s__")] += 1
                #################################
                #OPTION 2  -  confidence level #2
                #################################
                # if there is no entry in the provided taxonomy file, it looks for a  100% match, at >75% of the alignment length to anything
                # in the provided taxonomy file and assigns that taxonomy
                if success == "0":
                    rel_abund = 100*(float(line.split("\t")[1])/total)
                    blast_success = "0"
                    seq = ref_seqs_dict[line.split("\t")[0]]
                    sequence = seq.seq
                    with open("temp_seq_file.fasta", 'w') as temp_file:
                        temp_file.write(">temp_seq\n" + str(sequence))
                        temp_file.close
                    sub_string = "blastn -query temp_seq_file.fasta -db ~/Dropbox/QIIME/Latha/Graeme/utax_sequence_file.fasta -task blastn -dust no -outfmt \"7 qseqid sseqid evalue bitscore pident qcovhsp\" -perc_identity 100 -qcov_hsp_perc 75 -max_target_seqs 50 -num_threads 4 -out temp_blast_out.xml"
                    subprocess.call(sub_string, shell=True)
                    os.remove("temp_seq_file.fasta")
                    blast_OTU_success = "0"
                    # if the blast search has found anything:
                    if os.path.isfile("temp_blast_out.xml"):
                        if(int(file_len("temp_blast_out.xml")) > 5):
                            blast_OTU_success = "1"
                            with open("temp_blast_out.xml", 'r') as blast_output:
                                for outputline in blast_output:
                                    if outputline.startswith("temp_seq"):
                                        entrez_ID = outputline.split("\t")[1]
                                        perc_sim = outputline.split("\t")[4]
                                        hsp_perc = outputline.split("\t")[5]
                                sub_string2 = "grep \"" + entrez_ID + "\" " + args.input1 + " > temp_grep_output.txt"
                                subprocess.call(sub_string2, shell=True)
                            if os.path.isfile("temp_grep_output.txt"):
                                with open("temp_grep_output.txt", 'r') as grep:
                                    for result in grep:
                                        result = result.replace("\n", "")
                                        output_string = '\t'.join([result.rstrip("\n"),
                                                            "Percent_similarity:",
                                                            perc_sim,
                                                            "Percent_alignment_coverage",
                                                            hsp_perc.rstrip("\n")])
                                os.remove("temp_blast_out.xml")

                        # write out all the data from stage II of the process
                    if blast_OTU_success == "1":
                        output_file.write('\t'.join([output_string.split("\t")[0],
                                                    line.split("\t")[1].rstrip("\n"),
                                                    str(rel_abund),
                                                    output_string.split("\t")[1],
                                                    output_string.split("\t")[2],
                                                    output_string.split("\t")[3],
                                                    output_string.split("\t")[4],
                                                    output_string.split("\t")[5],
                                                    str(sequence)]) + "\n")
                #
                # OPTION 3  -  anything which still does not have a taxonomy:
                # the sequence is written out into a fasta file and blasted again the entire ncbi nt database
                # blast paramaters are for 95% sequence similarity for 95% of
                    if blast_OTU_success == "0":
                        list_of_data_to_blast.append('\t'.join([line.split("\t")[0],
                                                                line.split("\t")[1].rstrip("\n"),
                                                                str(rel_abund),
                                                                str(sequence)
                                                                ]))
                        total_lengths = total_lengths + int(len(str(sequence)))
                        list_of_fragment_sizes.append(int(len(str(sequence))))
                        blast_input.write(">" + line.split("\t")[0] + "\n" + str(sequence) + "\n")

# calculate average length of sequence to blast
average = int(total_lengths) / int(len(list_of_fragment_sizes))
coverage = float(average) / float(args.size_of_ref_fragment)
perc_coverage = float(coverage) * 100

## do the main round of blast searches on everything from step 3.
# this take a file generated previously "temo_blast_input.fasta", performs all the blast queries and writes out "temp_blast_out.xml"

sub_string = "blastn -query temp_blast_input.fasta -db ~/ncbi_nt/all_nt_merged -task blastn -dust no -outfmt \"7 qseqid sseqid evalue bitscore pident qcovhsp sscinames staxids\" -perc_identity 95 -qcov_hsp_perc  " + str(perc_coverage) + " -max_target_seqs 50 -num_threads 8 -out temp_blast_out2.xml"
subprocess.call(sub_string, shell=True)

# get the top result from each blast search
# split data rows into those which have a blast result and those which do not
blast_results = []
if os.path.isfile("temp_blast_out2.xml"):
    with open("temp_blast_out2.xml") as blast_output:
        for line in blast_output:
            new_result = 0
            if line.rstrip("\n") == "# BLASTN 2.2.31+":
                new_result = 1
                for i in range(0, 3):
                    line = blast_output.next()
                if line.rstrip("\n") == "# 0 hits found":
                    blast_results.append("No BLAST result found.")
                else:
                    line = blast_output.next()
                    line = blast_output.next()
                        # parse the data
                    ncbi_id = line.split("\t")[1]
                    blast_results.append('\t'.join([line.split("\t")[6],
                                                ncbi_id.split("|")[1],
                                                line.split("\t")[2],
                                                line.split("\t")[3],
                                                line.split("\t")[4],
                                                line.split("\t")[5],
                                                line.split("\t")[0]]) + "\n")
    os.remove("temp_blast_out2.xml")

data_present = []
empty_data = []
get_taxonomy = []
with open("temp_taxo_file.txt", 'w') as taxo:
    for x, y in zip(list_of_data_to_blast, blast_results):
        if y == "No BLAST result found.":
            empty_data.append('\t'.join([x.split("\t")[0],
                                        x.split("\t")[1],
                                        x.split("\t")[2],
                                        y.rstrip("\n"),
                                        x.split("\t")[3]
                                        ]))
        else:
            data_present.append('\t'.join([x.split("\t")[0],
                                           x.split("\t")[1],
                                           x.split("\t")[2],
                                           y,
                                           x.split("\t")[3]
                                           ]))
            taxo.write(y.split("\t")[1].rstrip("\n") + "\t999999999\n")

### write taxonomy IDs out into temporary file
# required for the gb_taxonomy_tools packages to run correctly(https://github.com/spond/gb_taxonomy_tools)
## these can now be checked against the ncbi taxonomy database to retrieve the taxo info

command_string = "gid-taxid temp_taxo_file.txt " + args.taxidnucl + " > gid-taxid_output.txt"
subprocess.call(command_string, shell=True)
os.remove("temp_taxo_file.txt")
# the ids generated above can then be checked against the main database to get the long-format taxonomy info
command_string = "cat gid-taxid_output.txt | taxonomy-reader " + args.names, \
                        + " " + args.nodes + " > taxonomy_output.txt"
subprocess.call(command_string, shell=True)

broken_taxonomy_output = []
taxonomy_info = []
if os.path.isfile("taxonomy_output.txt"):
    with open("taxonomy_output.txt", 'r') as taxo_file:
        for item in taxo_file:
            # see if there are any strings of 'n's in the taxonomy data
            # remove any that re found
            # this really destroys the output formatting if not done completely
            # if the taxonomy info is just a row full on N (some of them are)
            n_count = 0
            for x in item.split("\t"):
                if x == "n":
                    n_count = n_count + 1
            if n_count == 19:
                taxonomy_info.append('\t'.join([item.split("\t")[0],
                                    item.split("\t")[1],
                                    item.split("\t")[2],
                                    "No taxonomy information found"
                                    ]))
            else:
                item = item.replace("\tn\tn\tn", "\tn\tn")
                item = item.replace("\tn\tn", "\tn")
                taxonomy_info.append(item)

    os.remove("taxonomy_output.txt")
    output = []
    for item, row in zip(data_present, taxonomy_info):
        row = row.replace("\n", "")
        if row.split("\t")[3] == "No taxonomy information found" or len(row.split("\t")) < 16:
            broken_taxonomy_output.append('\t'.join([item.split("\t")[0],
                                    item.split("\t")[1],
                                    item.split("\t")[2],
                                    "No taxonomy information found in Blast , \
                                        Database. Blast result gave species as ", \
                                         + item.split("\t")[3],
                                    "Percent_similarity:",
                                    item.split("\t")[7].rstrip("\n"),
                                    "Percent_alignment_coverage:",
                                    item.split("\t")[8].rstrip("\n"),
                                    item.split("\t")[10].rstrip("\n"),
                                    "\n"
                                    ]))
        else:
            output.append('\t'.join([item.split("\t")[0],
                                    item.split("\t")[1],
                                    item.split("\t")[2],
                                    "k__" + row.split("\t")[5] + "; p__", \
                                        + row.split("\t")[7] + "; c__", \
                                         + row.split("\t")[9] + "; o__", \
                                          + row.split("\t")[11] + "; f__", \
                                           + row.split("\t")[13] + "; g__", \
                                            + row.split("\t")[15],
                                    "Percent_similarity",
                                    item.split("\t")[7],
                                    "Percent_alignment_coverage",
                                    item.split("\t")[8].rstrip("\n"),
                                    item.split("\t")[10].rstrip("\n"),
                                    "\n"
                                    ]))
with open(args.output, 'a') as output_file:
    for x in output:
        output_file.write(x)
    for x in broken_taxonomy_output:
        output_file.write(x)
    for x in empty_data:
        output_file.write(x + "\n")

KINGDOM_COUNTER2 = {}
PHYLA_COUNTER2 = {}
CLASS_COUNTER2 = {}
ORDER_COUNTER2 = {}
FAMILY_COUNTER2 = {}
GENUS_COUNTER2 = {}
SPECIES_COUNTER2 = {}

# get length of file
file_length =  file_len(args.output)

with open(args.output, 'r') as input_file:
    next(input_file)
    count = 1
    for line in input_file:
        print(line)
        if int(count) < int(file_length - 1):
            taxa = line.split("\t")[3]
            matchObj1 = re.findall( r'No BLAST result found', taxa, re.M|re.I)
            matchObj2 = re.findall( r'No taxonomy information found in Blast Database', taxa, re.M|re.I)
            if (matchObj1 or matchObj2):
                continue
            else:
                if not taxa.split("; ")[0].lstrip("k__") in KINGDOM_COUNTER2:
                    KINGDOM_COUNTER2[taxa.split("; ")[0].lstrip("k__")] = 1
                else:
                    KINGDOM_COUNTER2[taxa.split("; ")[0].lstrip("k__")] += 1
                if not taxa.split("; ")[1].lstrip("p__") in PHYLA_COUNTER2:
                    PHYLA_COUNTER2[taxa.split("; ")[1].lstrip("p__")] = 1
                else:
                    PHYLA_COUNTER2[taxa.split("; ")[1].lstrip("p__")] += 1
                if not taxa.split("; ")[2].lstrip("c__") in CLASS_COUNTER2:
                    CLASS_COUNTER2[taxa.split("; ")[2].lstrip("c__")] = 1
                else:
                    CLASS_COUNTER2[taxa.split("; ")[2].lstrip("c__")] += 1
                if not taxa.split("; ")[3].lstrip("o__") in ORDER_COUNTER2:
                    ORDER_COUNTER2[taxa.split("; ")[3].lstrip("o__")] = 1
                else:
                    ORDER_COUNTER2[taxa.split("; ")[3].lstrip("o__")] += 1
                if not taxa.split("; ")[4].lstrip("f__") in FAMILY_COUNTER2:
                    FAMILY_COUNTER2[taxa.split("; ")[4].lstrip("f__")] = 1
                else:
                    FAMILY_COUNTER2[taxa.split("; ")[4].lstrip("f__")] += 1
                if not taxa.split("; ")[5].lstrip("g__") in GENUS_COUNTER2:
                    GENUS_COUNTER2[taxa.split("; ")[5].lstrip("g__")] = 1
                else:
                    GENUS_COUNTER2[taxa.split("; ")[5].lstrip("g__")] += 1
                #if not taxa.split("; ")[6].lstrip("s__") in SPECIES_COUNTER2:
                #    SPECIES_COUNTER2[taxa.split("; ")[6].lstrip("s__")] = 1
                #else:
                #    SPECIES_COUNTER2[taxa.split("; ")[6].lstrip("s__")] += 1
            count = count + 1

with open(args.summary_output, 'w') as summary_out:
    summary_out.write("Top section of taxa info is derived only from the data that had entires in the utax taxonomy info. ie. from the closed ref OTU picking\n")
    summary_out.write("it does not include anything from the de-novo picking, although some (many) of those OTUs likely fall into the taxa listed here\n")
    summary_out.write("\nTAXA\tCOUNT\n")
    summary_out.write("\nKINGDOM\n")
    for x in KINGDOM_COUNTER:
        summary_out.write(x + "\t" + str(KINGDOM_COUNTER[x]) + "\n")
    summary_out.write("\nPHYLA\n")
    for x in PHYLA_COUNTER:
        summary_out.write(x + "\t" + str(PHYLA_COUNTER[x]) + "\n")
    summary_out.write("\nCLASS\n")
    for x in CLASS_COUNTER:
        summary_out.write(x + "\t" + str(CLASS_COUNTER[x]) + "\n")
    summary_out.write("\nORDER\n")
    for x in ORDER_COUNTER:
        summary_out.write(x + "\t" + str(ORDER_COUNTER[x]) + "\n")
    summary_out.write("\nFAMILY\n")
    for x in FAMILY_COUNTER:
        summary_out.write(x + "\t" + str(FAMILY_COUNTER[x]) + "\n")
    summary_out.write("\nGENUS\n")
    for x in GENUS_COUNTER:
        summary_out.write(x + "\t" + str(GENUS_COUNTER[x]) + "\n")
    #summary_out.write("\nSPECIES\n")
    #for x in SPECIES_COUNTER:
    #    summary_out.write(x + "\t" + str(SPECIES_COUNTER[x]) + "\n")

    summary_out.write("\n\n\nFollowing section contains the taxa info gathered from all different processing steps\n")
    summary_out.write("This includes the de novo OTU picking and the info gained from blasting the ncbi\n")
    summary_out.write("\nTAXA\tCOUNT\n")
    summary_out.write("\nKINGDOM\n")
    for x in KINGDOM_COUNTER2:
        summary_out.write(x + "\t" + str(KINGDOM_COUNTER2[x]) + "\n")
    summary_out.write("\nPHYLA\n")
    for x in PHYLA_COUNTER2:
        summary_out.write(x + "\t" + str(PHYLA_COUNTER2[x]) + "\n")
    summary_out.write("\nCLASS\n")
    for x in CLASS_COUNTER2:
        summary_out.write(x + "\t" + str(CLASS_COUNTER2[x]) + "\n")
    summary_out.write("\nORDER\n")
    for x in ORDER_COUNTER2:
        summary_out.write(x + "\t" + str(ORDER_COUNTER2[x]) + "\n")
    summary_out.write("\nFAMILY\n")
    for x in FAMILY_COUNTER2:
        summary_out.write(x + "\t" + str(FAMILY_COUNTER2[x]) + "\n")
    summary_out.write("\nGENUS\n")
    for x in GENUS_COUNTER2:
        summary_out.write(x + "\t" + str(GENUS_COUNTER2[x]) + "\n")
    #summary_out.write("\nSPECIES\n")
    #for x in SPECIES_COUNTER2:
    #    summary_out.write(x + "\t" + str(SPECIES_COUNTER2[x]) + "\n")
summary_out.close()
# tidy up some files
if os.path.isfile("temp_blast_input.fasta"):
    os.remove("temp_blast_input.fasta")
if os.path.isfile("temp_grep_output.txt"):
    os.remove("temp_grep_output.txt")
if os.path.isfile("gid-taxmid_output.txt"):
    os.remove("gid-taxid_output.txt")

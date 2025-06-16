# !/usr/bin/env python

#pass in list of proportions, identify proportion that is greatest, and return index of that proportion in list
def get_largest_proportion_idx(proportions_list):

    curr_max = -1
    curr_max_idx = -1

    for idx in range(0, len(proportions_list)):

        proportion = float(proportions_list[idx])

        if proportion > curr_max:

            curr_max = proportion
            curr_max_idx = idx

    return curr_max_idx

#get alleles present at position and put into list
def get_alleles(ref, alt):

    alt_alleles = alt.split(",")

    allele_list = [ref]
    allele_list += alt_alleles

    return allele_list

import sys, os, re

output_path = sys.argv[1]
sample_name = sys.argv[2]

output_file = open(output_path, "w")

output_file.write(">" + sample_name + "\n")

depth_search_pattern = "DP=(\d+)"
proportion_search_pattern = "QS=([^;]+)"

#go through each line of mpileup file
for line in sys.stdin:

    #ignore lines beginning with "#"
    if line[0] != "#":

        #split tab-delimited line
        split_line = line.strip().split("\t")

        #get chromosome, nucleotide position, reference allele, alternative allele, and info string
        chrom = split_line[0]
        pos = split_line[1]
        ref = split_line[3]
        alt = split_line[4]
        info_line = split_line[7]

        #skip lines that are marked as INDEL
        if info_line[:5] != "INDEL": 
            
            #get depth of reads at position
            depth_search = re.search(depth_search_pattern, info_line)
            depth = int(depth_search.group(1))

            #check if line has depth greater than or equal to 20
            if depth >= 20:

                #get the proportion of how frequently each allele appears, based on number of reads with that allele
                proportion_search = re.search(proportion_search_pattern, info_line)
                proportions = proportion_search.group(1)
                #put proportions into list
                proportions_list = proportions.split(",")

                #get list of alleles found at this position
                allele_list = get_alleles(ref, alt)

                #check that number of alleles is same as number of proportions
                if len(proportions_list) != len(allele_list):

                    print("something went wrong...", flush=True)
                    break

                #get index in list of greatest proportion
                dominant_idx = get_largest_proportion_idx(proportions_list)

                #use index to get allele present in greatest proportion
                dominant_allele = allele_list[dominant_idx]

                #write dominant allele to consensus sequence
                output_file.write(dominant_allele)

            #if line has depth less than 20
            else:

                #add to consensus the same number of N's as length of reference 
                for i in range(0, len(ref)):

                    output_file.write("N")

output_file.write("\n")

output_file.close()
#!/usr/bin/env python

import sys, os, math, time
from statistics import median

#determine window of read by doing floor of (read position / 100)
def determine_chromosome_window(split_read):

    read_pos = split_read[3]

    window = math.floor(int(read_pos)/100)

    return window

#update flag dictionary to count one more deletion
def increment_deletions(flag_dict):

    chromosome_window = determine_chromosome_window(split_read)

    if chromosome_window not in flag_dict:
        flag_dict[chromosome_window] = [0,1,0]
    else:
        flag_dict[chromosome_window][1] = flag_dict[chromosome_window][1] + 1
            
    return flag_dict

#update flag dictionary to count one more duplication
def increment_duplications(flag_dict):

    chromosome_window = determine_chromosome_window(split_read)
    
    if chromosome_window not in flag_dict:
        flag_dict[chromosome_window] = [0,0,1]
    else:
        flag_dict[chromosome_window][2] = flag_dict[chromosome_window][2] + 1

    return flag_dict

working_dir = sys.argv[1]
sample_name = sys.argv[2]
ena_name = sys.argv[3]
country = sys.argv[4]
sam_file_path = sys.argv[5]
output_path = sys.argv[6]
del_dup_bam_dir_path = sys.argv[7]

print("piped to flag_search_by_window.py")

os.system("mkdir -p " + working_dir)

### for extracting flagged reads ###

flagged_reads_sam_path = working_dir + sample_name + "_flagged_reads.sam"
os.system("samtools view -H " + sam_file_path + " -o " + flagged_reads_sam_path)
flagged_reads_sam_file = open(flagged_reads_sam_path, "a")

#####################################

bam_file = sys.stdin

#key is chromosome_window, value is list of int with length 3
#list = [#all, #deletion, #duplication]
flag_dict = {}

mapped_read_count = 0
del_flag_count = 0
dupl_flag_count = 0

bam_parse_start_time = time.time()
#iterate through reads
for read in bam_file:

    split_read = read.split("\t")
    chromosome = split_read[2]
    read_pos = split_read[3]

    #get flag and calculate window
    flag = int(split_read[1])

    if flag == 77 or flag == 141 or "Transfer" in chromosome or "MIT" in chromosome or "final" in chromosome:
        continue

    mapped_read_count += 1

    chromosome_window = determine_chromosome_window(split_read)

    #add to counter for all reads
    if chromosome_window not in flag_dict:
        flag_dict[chromosome_window] = [1,0,0]
    else:
        flag_dict[chromosome_window][0] = flag_dict[chromosome_window][0] + 1

    #check for flags for deletion and duplications
    if flag == 161 or flag == 81 or flag == 145 or flag == 97:

        #get index length
        index_length = int(split_read[8])

        #deletion
        if flag == 161 and index_length >= 1000:

            flag_dict = increment_deletions(flag_dict)
            del_flag_count += 1                        
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################

        elif flag == 161 and index_length <=-1000:

            flag_dict = increment_duplications(flag_dict)
            dupl_flag_count += 1
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################

        elif flag == 81 and index_length >= 1000:

            # DO NOT INCLUDE TO COUNT PAIRED ENDS
            #flag_dict = increment_duplications(flag_dict)
            #dupl_flag_count += 1
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################

        elif flag == 81 and index_length <= -1000:

            # DO NOT INCLUDE TO COUNT PAIRED ENDS
            #flag_dict = increment_deletions(flag_dict)
            #del_flag_count += 1                        
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################

        elif flag == 97 and index_length >= 1000:

            # DO NOT INCLUDE TO COUNT PAIRED ENDS
            #flag_dict = increment_deletions(flag_dict)
            #del_flag_count += 1                        
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################

        elif flag == 97 and index_length <= -1000:

            # DO NOT INCLUDE TO COUNT PAIRED ENDS
            #flag_dict = increment_duplications(flag_dict)
            #dupl_flag_count += 1
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################

        elif flag == 145 and index_length >= 1000:

            flag_dict = increment_duplications(flag_dict)
            dupl_flag_count += 1
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################

        elif flag == 145 and index_length <= -1000:
            
            flag_dict = increment_deletions(flag_dict)
            del_flag_count += 1                        
            ### for extracting flagged reads ###
            flagged_reads_sam_file.write(read)                        
            ####################################                  

### for extracting flagged reads ###
flagged_reads_sam_file.close()
####################################

bam_parse_end_time = time.time()
print("time to parse bam file: " + str(bam_parse_end_time-bam_parse_start_time), flush=True)

output_file = open(output_path, "a")

print("Flag_dict: " + str(flag_dict), flush=True)

summary_start_time = time.time()
#initialize coverage values
cov_min = 10000000000
cov_max = 0
cov_list = []

#initailize deletion values
del_min = 10000000000
del_max = 0
del_list = []

#initalize duplication values
dup_min = 10000000000
dup_max = 0
dup_list = []

del_window_list = []
dup_window_list = []

#iterate through flag dictionary and find min/median/max for coverage, deletion, duplication
for window in flag_dict:

    read_count = flag_dict[window][0]
    del_count = flag_dict[window][1]
    dup_count = flag_dict[window][2]

    #calculate coverage - forgot to chnge before running samples
    cov = (read_count * 200)/100

    #update list, use to calculate median later
    cov_list.append(cov)
    del_list.append(del_count)
    dup_list.append(dup_count)

    #update minimums
    if cov < cov_min:
        cov_min = cov
    if del_count < del_min:
        del_min = del_count
    if dup_count < dup_min:
        dup_min = dup_count

    #update maximums
    if cov > cov_max:
        cov_max = cov
    if del_count > del_max:
        del_max = del_count
    if dup_count > dup_max:
        dup_max = dup_count

    if del_count >= 5:
        del_window_list.append(str(window) + ":" + str(del_count))
    
    if dup_count >= 5:
        dup_window_list.append(str(window) + ":" + str(dup_count))


#calculate medians
cov_median = median(cov_list)
del_median = median(del_list)
dup_median = median(dup_list)

#write to output
dup_window_string = ""
del_window_string = ""

if len(del_window_list) > 0:
    del_window_string = ",".join(del_window_list)
else:
    del_window_string = "NA"

if len(dup_window_list) > 0:
    dup_window_string = ",".join(dup_window_list)
else:
    dup_window_string = "NA"

output_file.write(sample_name + "\t" + str(mapped_read_count) + "\t" + str(del_flag_count) + "\t" + str(dupl_flag_count) + "\t" + str(cov_min) + "\t" + str(cov_median) + "\t" + str(cov_max) + "\t" + str(del_min) + "\t" + str(del_median) + "\t" + \
                  str(del_max) + "\t" + str(dup_min) + "\t" + str(dup_median) + "\t" + str(dup_max) + "\t" +  del_window_string + "\t" + dup_window_string + "\n")

output_file.close()

if del_max >= 5 or dup_max >= 5:

    sort_idx_start_time = time.time()

    bam_file_path = del_dup_bam_dir_path + sample_name + ".sorted.bam"
    bai_file_path = del_dup_bam_dir_path + sample_name + ".sorted.bai"
    os.system("samtools sort " + sam_file_path + " -o " + bam_file_path)
    os.system("samtools index " + bam_file_path + " " + bai_file_path)

    ### for extracting flagged reads ###

    flagged_bam_file_path = working_dir + sample_name + "_flagged_reads.sorted.bam"
    flagged_bai_file_path = working_dir + sample_name + "_flagged_reads.sorted.bai"
    os.system("samtools sort " + flagged_reads_sam_path + " -o " + flagged_bam_file_path)
    os.system("samtools index " + flagged_bam_file_path + " " + flagged_bai_file_path)
    os.system("rm " + flagged_reads_sam_path)

    #####################################

    sort_idx_end_time = time.time()
    print("time to sort and index: " + str(sort_idx_end_time-sort_idx_start_time), flush=True)

summary_end_time = time.time()
print("time to calculate statistics: " + str(summary_end_time-summary_start_time), flush=True)

print(sample_name + " run done\n")


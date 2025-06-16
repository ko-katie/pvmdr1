#!/usr/bin/env python

import sys, os, re, time

start_time = time.time()
setup_start_time = time.time()

#read in arguments
working_dir = sys.argv[1]
samples_file_path = sys.argv[2]

#make subdirectories
mapping_dir = working_dir + "mapping_dir/"
flagged_reads_dir = working_dir + "flagged_reads_dir/"
del_dup_bam_files_dir = working_dir + "del_dup_bam_dir/"

os.system("mkdir -p " + working_dir + " " + mapping_dir + " " + del_dup_bam_files_dir)

#open output file and samples file
output_file_path = working_dir + "output.txt"
output_file = open(output_file_path, "w")

#add output file header
output_file.write("Sample\tReads Total\tDel Flags Total\tDupl Flags Total\tCov Min\tCov Median\tCov Max\tDel Min\tDel Median\tDel Max\tDup Min\tDup Median\tDup Max\tDel List\tDup List\n")
#output_file.write("Sample\tCov Min\tCov Median\tCov Max\tDel Min\tDel Median\tDel Max\tDup Min\tDup Median\tDup Max\n")
output_file.close()

#open bad_paths_file for paths that do not run successfully
bad_ena_file_path = working_dir + "bad_ena.txt"
bad_ena_file = open(bad_ena_file_path, "w")

print("started reading through samples", flush=True)

samples_file = open(samples_file_path, "r")
samples_file.readline()

setup_end_time = time.time()
print("time for setup: " + str(setup_end_time-setup_start_time), flush=True)

#iterate through each line of sample files to get r1 and r2 paths
for line in samples_file:

    sample_start_time = time.time()

    #split line by tab
    split_line = line.split("\t")
    sample_name = split_line[0]
    country = split_line[4].replace(" ", "")
    sample_ena_list = split_line[8].strip("\"").replace(" ", "")
    qc_pass = split_line[12]

    #skip samples with no accessions
    if len(sample_ena_list) == 0:
        print(str(sample_name) + " has no accessions")
        continue    

    if qc_pass == "TRUE":

        print("Sample name: " + str(sample_name))

        #initialize sam path for mapping
        sam_path = mapping_dir + sample_ena_list + ".sam"

        try:

            #mdr1 gene is from PvP01_10_v2:476,677..483,934, interested in 60kb region surrounding
            #make 100bp windows from 450000-510000
            mapping_start_time = time.time()
            print("perform mapping", flush=True)
            os.system("/usr/local/packages/hisat2/hisat2 -x /local/projects-t3/SerreDLab-3/kko/Pv_P01_Index/PvivaxP01 -f \
                --sra-acc " + sample_ena_list + " --max-intronlen 20 -S " + sam_path + " -p 12")

            mapping_end_time = time.time()
            print("time to map: " + str(mapping_end_time - mapping_start_time), flush=True)

        except:

            bad_ena_file.write(sample_ena_list + "\t")
            print(sample_ena_list + " could not be retrieved", flush=True)
            continue

        #pass sam file to flag_search_by_window.py script to count flags
        
        os.system("samtools view " + sam_path + " | python3 /local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/flag_search_by_window_final_sra.py " + flagged_reads_dir \
                + " " + sample_name + " " + sample_ena_list + " " + country + " " + sam_path + " " + output_file_path + " " + del_dup_bam_files_dir)

        #remove sam file
        os.system("rm " + sam_path)

        sample_end_time = time.time()
        print("Sample run time: " + str(sample_end_time-sample_start_time))

    else:
        print(str(sample_name) + " was not run")
        

#close files
samples_file.close()
bad_ena_file.close()

end_time = time.time()
print("Elapsed time: " + str(end_time-start_time))
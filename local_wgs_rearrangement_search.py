#!/usr/bin/env python

import sys, os, re, time

start_time = time.time()
setup_start_time = time.time()

#read in arguments
working_dir = sys.argv[1]
r1_samples_file_path = sys.argv[2]
r2_samples_file_path = sys.argv[3]

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
output_file.close()

#open bad_paths_file for paths that do not run successfully
bad_ena_file_path = working_dir + "bad_ena.txt"
bad_ena_file = open(bad_ena_file_path, "w")

print("started reading through samples", flush=True)

r1_samples_file = open(r1_samples_file_path, "r")
r2_samples_file = open(r2_samples_file_path, "r")

setup_end_time = time.time()
print("time for setup: " + str(setup_end_time-setup_start_time), flush=True)

sample_name_pattern = "UMB/(.+)/ILL"

#iterate through each line of sample files to get r1 and r2 paths
for r1_line in r1_samples_file:

    sample_start_time = time.time()

    r1_line_stripped = r1_line.strip()
    r2_line_stripped = r2_samples_file.readline().strip()

    if os.path.isfile(r1_line_stripped) and os.path.isfile(r2_line_stripped):

        regex_search = re.search(sample_name_pattern, r1_line)
        sample_name = regex_search.group(1)

        print("Sample name: " + str(sample_name))

        #initialize sam path for mapping
        sam_path = mapping_dir + sample_name + ".sam"

        try:

            #mdr1 gene is from PvP01_10_v2:476,677..483,934, interested in 20kb region surrounding
            #make 1kb windows from 470000-490000
            #change to hisat with --disable-splicer
            mapping_start_time = time.time()
            print("perform mapping", flush=True)
            os.system("/usr/local/packages/hisat2/hisat2 -x /local/projects-t3/SerreDLab-3/kko/Pv_P01_Index/PvivaxP01 -q \
                -1 " + r1_line_stripped +  " -2 " + r2_line_stripped + " --max-intronlen 20 -S " + sam_path + " -p 24")

            mapping_end_time = time.time()
            print("time to map: " + str(mapping_end_time - mapping_start_time), flush=True)

        except:

            bad_ena_file.write(sample_name + "\t")
            print(sample_name + " could not be retrieved", flush=True)
            continue

        #pass sam file to flag_search_by_window.py script to count flags
        os.system("samtools view " + sam_path + " | python3 /local/projects-t3/SerreDLab-3/kko/mdr_del_search_dna_whole_genome/flag_search_by_window_final_local.py " + flagged_reads_dir \
                + " " + sample_name + " " + sam_path + " " + output_file_path + " " + del_dup_bam_files_dir)

        #remove sam file
        os.system("rm " + sam_path)

        sample_end_time = time.time()
        print("Sample run time: " + str(sample_end_time-sample_start_time))

    else:
        print(str(sample_name) + " was not run")
        

#close files
r1_samples_file.close()
r2_samples_file.close()
bad_ena_file.close()

end_time = time.time()
print("Elapsed time: " + str(end_time-start_time))
#!/usr/bin/python
# Create a key for paired end read fastq files
#
# Pulls the fastq files present within a parent directory
# 'fastq_dir' specifies the parent directory.
# Within fastq_dir, each directory must correspond to a sample. The
# name of each directory is assumed to be the sample name
#
# The script recursively lists files in each sample folder, and grabs
# fastq files based on the suffix
#

import os
import sys
import gzip

def isfastq(x):
    if (x[-9:] == ".fastq.gz") or \
    (x[-6:] == ".fastq") or \
    (x[-3:] == ".fq") or \
    (x[-6:] == ".fq.gz"):
        return True
    else:
        return False

def isr1(x, r = ["R1_0", "R1001", "_1.fastq"]):
    for i in r:
        if i in x:
            return True
    return False

# def isr2(x, r = ["R2_0", "R2001", "_2.fastq"]):
#     for i in r:
#         if i in x:
#             return True
#     return False

# fastq_dir = "/u/project/pajukant/nikodm/sbc_sirna/data/fastq/"

# first argument is the path to the fastq directory
fastq_dir = sys.argv[1]
# second argument determines whether we want the raw or trimmed files
# if second argument is "raw", we want the raw files
raw = sys.argv[2] == 'raw'

samples = os.listdir(fastq_dir)

if (raw == False):
    sample_file_name = "../../data/sample_fastq.txt"
else:
    sample_file_name = "../../data/sample_fastq_raw.txt"
ofs = open(sample_file_name, 'w')
ofs.write("\t".join(["sample", "path"]) + "\n")

for sample in samples:
    fastq_dir_s = fastq_dir + sample
    for r, d, f in os.walk(fastq_dir_s):
        fastqs = [x for x in f if isfastq(x)]
        if (len(fastqs) == 0):
            continue
        f1files = sorted([x for x in fastqs if isr1(x)])
        if (raw == False):
            f1files = [x for x in f1files if "trimmed" in x]
        else:
            f1files = [x for x in f1files if "trimmed" not in x]
        #f2files = sorted([x for x in fastqs if isr2(x)])
        #if (len(f1files) != len(f2files)):
        #    print("Number of R1 fastq files doesn't match number of R2 fastq files")
        #    sys.exit(1)
        fastq1paths = ",".join([r + "/" + file for file in f1files])
        #fastq2paths = ",".join([r + "/" + file for file in f2files])

        fastqpaths = fastq1paths #+ ';' + fastq2paths
        lineout = "\t".join([sample, fastqpaths]) + "\n"
        ofs.write(lineout)
ofs.close()

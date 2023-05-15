#!/usr/bin/env python
#$ -S /u/local/apps/python/2.7.3/bin/python
#$ -cwd
#$ -j y
#$ -pe shared 12
#$ -l h_data=4000M,h_rt=3:00:00
#$ -v QQAPP=openmp
#$ -v LD_LIBRARY_PATH
#$ -v PYTHON_LIB
#$ -v PYTHON_DIR
#$ -m n
#$ -r n
#$ -o map.py.log

import os
import sys
from subprocess import call

p = sys.argv[1]

star="/u/project/pajukant/nikodm/bulk_practice_METSIM/bin/STAR-2.5.2b/bin/Linux_x86_64/STAR"
genomeDir = "../../data/genome_for_pass" + p

outdir = "../../data/aligned/pass" + p
if not os.path.exists(outdir):
        os.makedirs(outdir)

# Get Sample IDs
with open("../../data/sample_fastq.txt") as f:
  lines = f.readlines()

for line in lines[1:]:
  linel = line.replace(';', '\t').strip().split('\t')
  sampleID = linel[0]
  fastq_r1 = linel[1]
  # fastq_r2 = linel[2]
  read1paths = fastq_r1
  # read2paths = fastq_r2
  outpre = outdir + "/" + sampleID + "/" + sampleID + "_"
  if not os.path.exists(outpre):
    os.makedirs(outpre)
  call([star,
    "--runThreadN", "12",
    "--genomeDir", genomeDir,
    "--genomeLoad", "LoadAndKeep",
    "--readFilesIn", read1paths, #read2paths,
    "--readFilesCommand", "zcat",
    "--alignIntronMin", "20",
    "--alignIntronMax", "1000000",
    "--outFilterMultimapNmax", "1",
    "--chimSegmentMin", "15",
    "--outFilterMismatchNmax", "2",
    "--outFilterMismatchNoverLmax", "0.04",
    "--outSAMattributes", "All",
    "--outSAMtype", "BAM", "Unsorted",
    "--outFileNamePrefix", outpre])

call([star, "--genomeDir", genomeDir, "--genomeLoad", "Remove"])

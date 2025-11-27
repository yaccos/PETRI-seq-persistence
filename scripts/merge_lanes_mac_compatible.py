## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Tavazoie Lab, Columbia University
## Last updated March 2020 

#!/usr/bin/python

import sys
import glob
import re
import subprocess
import os

# combine all data from the different lanes of a nextseq run

fname_re = re.compile(r"(.*)_L00\d_R\d_001.fastq")


fastqfiles=glob.glob("*_001.fastq")
print(fastqfiles)

prefixes = set()

def get_prefix(infilename):
  # return the prfix of the read information in a given filename

  mymatch = fname_re.match(infilename)

  if mymatch is None:
    return None
 
  else:
    return mymatch.groups()[0]
  

for filename in fastqfiles:
  this_prefix = get_prefix(filename)

  if this_prefix is not None:
    prefixes.add(this_prefix)

for prefix in prefixes:
  print(f"Working on prefix {prefix}")
  these_files_r1 = glob.glob(f"{prefix}*_R1_*fastq")
  these_files_r2 = glob.glob(f"{prefix}*_R2_*fastq")
  these_files_r1.sort()
  these_files_r2.sort()
  file_str_1 = " ".join(these_files_r1)
  file_str_2 = " ".join(these_files_r2)
  print(f"cat {file_str_1} > {prefix}_R1_all_lanes.fastq")
  print(f"cat {file_str_2} > {prefix}_R2_all_lanes.fastq")
  subprocess.call(f"cat {file_str_1} > {prefix}_R1_all_lanes.fastq", shell=True)
  subprocess.call(f"cat {file_str_2} > {prefix}_R2_all_lanes.fastq", shell=True)
  for fname in these_files_r1:
    os.unlink(fname)

  for fname in these_files_r2:
    os.unlink(fname)


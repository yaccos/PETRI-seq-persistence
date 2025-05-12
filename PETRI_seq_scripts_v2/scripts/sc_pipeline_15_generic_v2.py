## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated January 2025

import time
import os
import sys
import rpy2.robjects as robjects

start = time.time()
script_dir = sys.argv[0].split(sys.argv[0].split('/')[-1])[0]
i = 1
sample = sys.argv[i][0:sys.argv[i].find('_S')]
if sys.argv[i].find('_S') == -1:
    raise ValueError('Must include _S after sample')
n_lanes = int(sys.argv[2])

print(f'Preprocessing {sample}')


robjects.globalenv['n_lanes'] = n_lanes
robjects.globalenv['script_dir'] = script_dir
robjects.globalenv['sample'] = sample

robjects.r.source(f'{script_dir}/demultiplexer.R')

end = time.time()
print(f"Time elapsed during Python pipeline: {end - start}")

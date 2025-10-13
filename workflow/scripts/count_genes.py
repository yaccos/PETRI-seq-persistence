import pysam
import pandas as pd
import sys
from dataclasses import dataclass
from typing import Dict, Tuple
import umi_tools
import multiprocessing
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

@dataclass(frozen=True)
class BarcodeGene:
    barcode: str
    gene: str
    contig: str


def get_gene(read: pysam.AlignedSegment):
    num_matches = int(read.get_tag("X0"))

    if read.has_tag('XT'):
        gene = read.get_tag("XT")
        if ('rRNA' not in gene) and (num_matches > 1):
            gene = 'ambiguous'
    else:
        gene = 'no_feature'
    return gene

def get_contig(read: pysam.AlignedSegment):
    proposed_contig = read.reference_name
    return proposed_contig

def get_alignment_status(read: pysam.AlignedSegment):
    if not read.is_mapped:
        contig = "unaligned"
        gene = "unaligned"
        return contig, gene
    contig = get_contig(read)
    gene = get_gene(read)
    return contig, gene

def prepare_barcode_table(file_path: str):
    barcode_table = pd.read_table(file_path).set_index("read").sort_index()
    barcode_table["celltag"] = barcode_table[["bc1", "bc2", "bc3"]].agg('_'.join, axis=1)
    return dict(zip(barcode_table.index, zip(barcode_table["UMI"], barcode_table["celltag"])))

def count_umis(umi_dict: Dict[bytes, int]):
    clusterer = umi_tools.UMIClusterer(cluster_method="directional")
    clustered_umis = clusterer(umi_dict, threshold=1)
    return len(clustered_umis)

def create_count_record(dict_pair: Tuple[BarcodeGene, Dict[bytes, int]]):
    barcode_info, umi_dict = dict_pair
    contig_gene = f"{barcode_info.contig}:{barcode_info.gene}"
    count = count_umis(umi_dict)
    return (barcode_info.barcode, contig_gene, count)


threshold = int(sys.argv[1])
sample = sys.argv[2]
n_cores = int(sys.argv[3])
chunk_size = int(sys.argv[4])
barcode_table_file = f"results/{sample}/{sample}_barcode_table.txt"
logging.info("Reading barcode table")
barcode_table = prepare_barcode_table(barcode_table_file)
bam_file_path = f"results/{sample}/{sample}_sorted.bam.featureCounts.bam"
bamfile = pysam.AlignmentFile(bam_file_path, "rb")
cell_UMI_count: Dict[BarcodeGene, Dict[bytes, int]] = {}
logging.info("Parsing reads")
iteration_gap = int(1e6)
read_count = 0
for read in bamfile.fetch(until_eof=True):
    read_name = str(read.query_name)
    if read_name not in barcode_table:
         # Read is filtered out
         continue
    read_count += 1
    contig, gene = get_alignment_status(read)
    if contig == 'ambiguous' or gene == "ambiguous":
         continue
    barcode_info = barcode_table[read_name]
    cell_barcode = str(barcode_info[1])
    # UMI-tools expect the UMI to be in the form of bytes
    cell_umi = bytes(barcode_info[0], encoding="ascii")
    read_key = BarcodeGene(cell_barcode, gene, contig)
    if read_key not in cell_UMI_count:
         cell_UMI_count[read_key] = {cell_umi: 1}
    elif cell_umi not in cell_UMI_count[read_key]:
         cell_UMI_count[read_key][cell_umi] = 1
    else:
         cell_UMI_count[read_key][cell_umi] += 1
    if read_count % iteration_gap == 0:
        logging.info(f"Parsed {read_count} reads")


logging.info(f"Deduplicating UMIs")
with multiprocessing.Pool(n_cores) as p:
    res_list = p.map(create_count_record, cell_UMI_count.items(), chunksize=chunk_size)

logging.info(f"Preparing results")

res_frame = pd.DataFrame.from_records(res_list, columns=["Cell Barcode","contig_gene","count"])

n_operons = res_frame["count"].sum()
mean_operons = res_frame["count"].mean()
max_operons = res_frame["count"].max()


(res_frame.loc[lambda df: df["count"] > threshold].
 pivot_table(index="Cell Barcode", columns="contig_gene", values="count", aggfunc='sum', fill_value=0).
 to_csv(f"results/{sample}/{sample}_gene_count_matrix.txt",sep='\t')
)

logging.info("DONE")
logging.info(f"Total number of reads: {read_count}")
logging.info(f"Total number of groups deduplicated: {len(cell_UMI_count)}")
logging.info(f"Total number of unique UMIs: {n_operons}")
logging.info(f"Mean number of unique UMIs per operon: {mean_operons}")
logging.info(f"Max number of unique UMIs per operon: {max_operons}")

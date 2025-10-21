import pysam
import pandas as pd
import sys
from dataclasses import dataclass
from typing import Dict, Tuple
import umi_tools
import multiprocessing
import logging
import sqlite3
from functools import lru_cache

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
summary_logger = logging.getLogger("count_genes.summary")
if not summary_logger.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(message)s"))
    summary_logger.addHandler(_handler)
summary_logger.propagate = False

@dataclass(frozen=True)
class BarcodeGene:
    barcode: str
    gene: str
    contig: str


def prepare_barcode_table(db_path: str):
    con = sqlite3.connect(db_path)
    cursor = con.cursor()
    query = "SELECT UMI, celltag FROM selected_barcodes WHERE read = ?"

    def lookup(read_id: str):
        row = cursor.execute(query, (read_id,)).fetchone()
        if row is None:
            return None
        return row  # (UMI, celltag)

    return con, lookup


def is_ambiguous(read: pysam.AlignedSegment) -> bool:
    if read.mapping_quality == 0:
        return True
    if read.has_tag("SA"):
        return True
    if read.has_tag("XG") and read.get_tag("XG") != 0:
        return True
    return False


def get_gene(read: pysam.AlignedSegment):
    if not read.has_tag("XT"):
        return "no_feature"
    gene: str = read.get_tag("XT")
    if "rRNA" not in gene and is_ambiguous(read):
        return "ambiguous"
    return gene

def get_contig(read: pysam.AlignedSegment):
    proposed_contig = read.reference_name
    return proposed_contig

def get_alignment_status(read: pysam.AlignedSegment):
    contig = get_contig(read)
    gene = get_gene(read)
    return contig, gene

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
barcode_table_file = f"results/{sample}/{sample}_selected_barcode_table.sqlite"
logging.info("Opening barcode database")
barcode_table = prepare_barcode_table(barcode_table_file)
barcode_con, fetch_barcode = prepare_barcode_table(barcode_table_file)
bam_file_path = f"results/{sample}/{sample}_sorted.bam.featureCounts.bam"
bamfile = pysam.AlignmentFile(bam_file_path, "rb")
cell_UMI_count: Dict[BarcodeGene, Dict[bytes, int]] = {}
logging.info("Parsing alignments")
iteration_gap = int(1e6)
alignment_count = 0
read_count = 0
selected_count = 0
aligned_count = 0
unambiguous_count = 0
feature_determined = 0

for read in bamfile.fetch(until_eof=True):
    alignment_count += 1
    if alignment_count % iteration_gap == 0:
        logging.info(f"Parsed {alignment_count} alignments")
    if read.is_secondary or read.is_supplementary:
        continue
    read_count += 1
    read_name = str(read.query_name)
    barcode_info = fetch_barcode(read_name)
    if barcode_info is None:
         # Read is filtered out
         continue
    selected_count += 1

    if not read.is_mapped:
        continue
    aligned_count += 1
    contig, gene = get_alignment_status(read)
    if gene == "ambiguous":
        continue
    unambiguous_count += 1
    if gene != "no_feature":
        feature_determined += 1
    
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


logging.info(f"Deduplicating UMIs")
with multiprocessing.Pool(n_cores) as p:
    res_list = p.imap_unordered(create_count_record, cell_UMI_count.items(),chunksize=chunk_size)
    res_list = list(res_list)

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

def round_percentage(total, part, ndigits=2):
    if total == 0:
        # In order to avoid zero division errors
        total = 1
    percentage = part / total * 100
    return round(percentage, ndigits)

logging.basicConfig(level=logging.INFO, format='%(message)s')

summary_logger.info("")
summary_logger.info(f"Total number of alignments processed: {alignment_count}")
summary_logger.info(f"of which {read_count} ({round_percentage(alignment_count, read_count)}%) alignments were primary alignments of a read:")
summary_logger.info(f"of which {selected_count} ({round_percentage(read_count, selected_count)}%) reads were selected based on barcode frequency")
summary_logger.info(f"of which {aligned_count} ({round_percentage(selected_count, aligned_count)}%) reads aligned to the reference genome")
summary_logger.info(f"of which {unambiguous_count} ({round_percentage(aligned_count, unambiguous_count)}%) reads aligned unambiguously to the reference genome")
summary_logger.info(f"of which {feature_determined} ({round_percentage(unambiguous_count, feature_determined)}%) reads did match a genomic feature")
summary_logger.info(f"Total number of groups deduplicated: {len(cell_UMI_count)}")
summary_logger.info(f"Total number of unique UMIs: {n_operons}")
summary_logger.info(f"Mean number of unique UMIs per operon: {mean_operons}")
summary_logger.info(f"Max number of unique UMIs per operon: {max_operons}")

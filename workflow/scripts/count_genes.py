import pysam
import pandas as pd
import sys
from dataclasses import dataclass
from typing import Dict, Tuple
import umi_tools

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
    if proposed_contig == "CP000255":
         pass
    if not read.has_tag("X0"):
        return proposed_contig
    num_matches = int(read.get_tag("X0"))
    if num_matches <= 1:
        return proposed_contig
    
    if not read.has_tag("XA"):
        return "ambiguous"
    
    edit_dist = read.get_tag("NM")
    match_list = read.get_tag("XA").split(';')

    for entry in match_list:
        if entry != '' and entry.split(',')[3] == edit_dist and entry.split(',')[0] not in proposed_contig:
                    return 'ambiguous'
        
    return proposed_contig

def get_alignment_status(read: pysam.AlignedSegment):
    if not read.is_mapped:
        contig = "unaligned"
        gene = "unaligned"
        return contig, gene
    contig = get_contig(read)
    gene = get_gene(read)
    return contig, gene

# sample = sys.argv[1]
sample = "random20000"
barcode_table_file = f"../results/{sample}/{sample}_barcode_table.txt"
barcode_table = pd.read_table(barcode_table_file).set_index("read").sort_index()
barcode_table["celltag"] = barcode_table[["bc1", "bc2", "bc3"]].agg('_'.join, axis=1)
bam_file_path = f"../results/{sample}/{sample}_sorted.bam.featureCounts.bam"
bamfile = pysam.AlignmentFile(bam_file_path, "rb")
cell_UMI_count: Dict[BarcodeGene, Dict[bytes, int]] = {}
for read in bamfile.fetch():
    read_name = str(read.query_name)
    
    if read_name not in barcode_table.index:
         # Read is filtered out
         continue
    contig, gene = get_alignment_status(read)
    if contig == 'ambiguous' or gene == "ambiguous":
         continue
    cell_barcode = str(barcode_table.loc[read_name, "celltag"])
    cell_umi = bytes(barcode_table.loc[read_name, "UMI"], "ascii")
    read_key = BarcodeGene(cell_barcode, gene, contig)
    if read_key not in cell_UMI_count:
         cell_UMI_count[read_key] = {cell_umi: 1}
    elif cell_umi not in cell_UMI_count[read_key]:
         cell_UMI_count[read_key][cell_umi] = 1
    else:
         cell_UMI_count[read_key][cell_umi] += 1

cell_UMI_summary: Dict[BarcodeGene, int] = {}
clusterer = umi_tools.UMIClusterer(cluster_method="directional")
for barcode_gene, umi_dict in cell_UMI_count.items():
    clustered_umis = clusterer(umi_dict, threshold=1)
    cell_UMI_summary[barcode_gene] = len(clustered_umis)

res_list = [(key.barcode, f"{key.contig}:{key.gene}", value) 
                 for key, value in cell_UMI_summary.items()] 
res_frame = (
     pd.DataFrame.from_records(res_list, columns=["Cell Barcode","contig_gene","count"]).
     pivot_table(index="Cell Barcode", columns="contig_gene", values="count", fill_value=0)
)

res_frame.to_csv(f"../results/{sample}/{sample}_gene_count_matrix.txt",sep='\t')

#sam processor for directional UMI grouping
import sys
import re

version = 11
threshold = int(sys.argv[1])
sample = sys.argv[2]


def get_tag(line, tag):
    match = re.search(f"(?<={tag}:\\S:)\\S+", line)
    return match.group()


def get_contig(line):

    proposed_contig = line.split('\t')[2] 

    if 'X0:i' not in line:
        return proposed_contig
    

    num_matches = int(get_tag(line, "X0"))
    if num_matches <= 1:
         return proposed_contig

    
    if 'XA:Z' not in line:
        return 'ambiguous'
    
    edit_dist = get_tag(line, "NM")
    match_list = get_tag(line, "XA").split(';')
    for entry in match_list:
        if entry != '' and entry.split(',')[3] == edit_dist and entry.split(',')[0] not in proposed_contig:
                    return 'ambiguous'
    return proposed_contig

def get_gene(line):
    num_matches = int(get_tag(line, "X0"))
    if 'XT:Z' in line:
        gene = get_tag(line, "XT")
        if gene == "NA":
           # This NA tag comes from add_cell_barcode.R which assigns NA values to non-existent tags
           gene = "no_feature"  
        elif ('rRNA' not in gene) and (num_matches > 1):
            gene = 'ambiguous'
    else:
        gene = 'no_feature'
    
    return gene
     

     

def get_alignment_status(line):

    if line.split('\t')[1] == '4':
            contig = 'unaligned'
            gene = 'unaligned'
            return contig, gene
    
    contig = get_contig(line)
    gene = get_gene(line)
    
    return contig, gene

UMIgene_to_count = {}

sam = open(f'results/{sample}_FC_directional_grouped_2/{sample}_group_FC.sam', 'r')

for line in sam:
        cell_barcode = get_tag(line, "CB")

        if cell_barcode not in UMIgene_to_count:
            UMIgene_to_count[cell_barcode] = {}
        
        barcode_count_dict = UMIgene_to_count[cell_barcode]

        UMI = get_tag(line, "BX")
        contig, gene = get_alignment_status(line)

        UMIgene = (UMI, contig, gene)

        if UMIgene in barcode_count_dict:
            barcode_count_dict[UMIgene] += 1
        else:
            barcode_count_dict[UMIgene] = 1

sam.close()

output = open(f'results/{sample}_v{version}_threshold_{threshold}_filtered_mapped_UMIs.txt', 'w')
output.write('Cell Barcode\tUMI\tcontig:gene\ttotal_reads\n')

for cell_barcode, count_dict in UMIgene_to_count.items():
    for UMIgene, count in count_dict.items():
        UMI = UMIgene[0]
        contig_gene = f'{UMIgene[1]}:{UMIgene[2]}'
        if count > threshold:
            output.write(f'{cell_barcode}\t{UMI}\t{contig_gene}\t{count}\n')

output.close()

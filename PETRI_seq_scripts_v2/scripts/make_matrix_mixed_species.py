import sys
import pandas as pd

sample = sys.argv[1]
table = pd.read_csv(sample + '_filtered_mapped_UMIs.txt',sep='\t',index_col=0)
gene_matrix = (table[['contig:gene','UMI']].
               groupby(['Cell Barcode','contig:gene']).
               count().
               unstack(level='Cell Barcode').
               fillna(0).
               transpose().
               droplevel(0).
               filter(regex="^(?!.*ambiguous).*$") # We filter away the columns *not* containing 'ambiguous' by using a negative lookahead
               )
gene_matrix.to_csv(sample + '_mixed_species_gene_matrix.txt',sep='\t')


import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
import scipy.stats as st
import random
random.seed(1)
from collections import Counter
from scipy.stats import norm
import sys
from scipy.stats.mstats import gmean

sample1 = 'SBC211'
sample2 = 'SBC212'
dirs1 = ['source_data/included/'] # indicate where crRNA count file is
dirs2 = ['source_data/included/'] # indicate where crRNA count file is
sample1_genomes = 'MG1655_pC2185_pBbS6Ccas9_NCTC8325_pWJ411_pWJ418_pWJ420' # table name includes all references used for alignment (E. coli genome, staph genome (from library generation), plasmids)
sample2_genomes = 'MG1655_pC2185_pBbS6Ccas9_NCTC8325_pWJ411_pWJ418_pWJ420' # table name includes all references used for alignment (E. coli genome, staph genome (from library generation), plasmids)
n_rep = 100 # number of times to re-sample the null distribution to get p-value; default set to 100000

def combine_LAG_tables_logfc(table1,table2,sample1,sample2):
    table1 = table1.add_suffix('.' + sample1)
    table1.columns = [table1.columns[0],sample1,table1.columns[2],table1.columns[3]]
    table2 = table2.add_suffix('.' + sample2)
    table2.columns = [table2.columns[0],sample2,table2.columns[2],table2.columns[3]]
    LAG_table = table1.merge(table2,left_index=True,right_index=True,how='inner')
    LAG_table = LAG_table[[str('chr_dir_pos.' + sample1),str('LAGs.' + sample1),str('LAG_orientations.' + sample1),sample1,sample2]]
    LAG_table.columns = ['chr_dir_pos','LAGs','orientation',sample1,sample2]
    samples = [sample1,sample2]
    for sample in samples:
        LAG_table['pseudocount.' + sample] = LAG_table[sample]
        LAG_table.loc[LAG_table['pseudocount.' + sample]==0,str('pseudocount.' + sample)]+=0.99
        LAG_table['freq.' + sample] = LAG_table['pseudocount.' + sample]/sum(LAG_table[sample])
    LAG_table['log2fc'] = np.log2(LAG_table['freq.' + sample2]/LAG_table['freq.' + sample1])
    return LAG_table
def filter_LAG(LAG_table,sample1,sample2,threshold):
    LAG_table_filtered = LAG_table.loc[(LAG_table[sample1]>=threshold)|(LAG_table[sample2]>=threshold)]
    return LAG_table_filtered
def expand_multi_LAGs(LAG_table):
    table = LAG_table.dropna()
    table['LAGs'] = table['LAGs'].str.split(';')
    table['orientation'] = [list(x) for x in table['orientation']]
    table['strand'] = [list(x) for x in table['strand']]
    table_expanded = table.explode('LAGs')
    table_expanded['orientation'] = table['orientation'].explode()
    table_expanded['strand'] = table['strand'].explode()
    return table_expanded

def add_sense(LAG_table):
    dirs = LAG_table['chr_dir_pos'].str.split('-',expand=True)[1]
    LAG_table.loc[((LAG_table['orientation']=='+')&(dirs=='0')) | ((LAG_table['orientation']=='-')&(dirs=='16')),'strand'] = 's'
    LAG_table.loc[((LAG_table['orientation']=='-')&(dirs=='0')) | ((LAG_table['orientation']=='+')&(dirs=='16')),'strand'] = 'a'
    multi_LAG = LAG_table.loc[LAG_table['strand'].isna(),:]
    sense_list = []
    for i in range(0,len(multi_LAG)):
        d = multi_LAG.iloc[i]['chr_dir_pos'].split('-')[1]
        o = list(multi_LAG.iloc[i,:]['orientation']) # orientation
        l = multi_LAG.iloc[i,:]['LAGs'].split(';') # LAG
        s = ''
        for char in o:
            if ((d=='0') & (char=='+')) | ((d=='16') & (char=='-')):
                s = s + 's'
            else:
                s = s + 'a'
        sense_list.append(s)
    LAG_table.loc[LAG_table['strand'].isna(),'strand'] = sense_list
    return LAG_table

counts = 0
for d in dirs2:
    file_name = d + sample2 + '_' + sample2_genomes + '_CDP_hit_LAG_ori_6.txt'
    if counts == 0:
        table = pd.read_csv(file_name,sep='\t',index_col=0)
        counts = 1
    else:
        table.iloc[:,0] = table.iloc[:,0] + pd.read_csv(file_name,sep='\t',index_col=0).iloc[:,0]
table['chr_dir_pos'] = table.index
table.index = range(0,len(table))
table2 = table.loc[:,pd.read_csv(file_name,sep='\t').columns]

counts = 0
for d in dirs1:
    file_name = d + sample1 + '_' + sample1_genomes + '_CDP_hit_LAG_ori_6.txt'
    if counts == 0:
        table = pd.read_csv(file_name,sep='\t',index_col=0)
        counts = 1
    else:
        table.iloc[:,0] = table.iloc[:,0] + pd.read_csv(file_name,sep='\t',index_col=0).iloc[:,0]
table['chr_dir_pos'] = table.index
table.index = range(0,len(table))
table1 = table.loc[:,pd.read_csv(file_name,sep='\t').columns]

threshold=10
all_LAGs = combine_LAG_tables_logfc(table1,table2,sample1,sample2)
all_LAGs_filtered = filter_LAG(all_LAGs,sample1,sample2,threshold)

NA_group_treated = all_LAGs_filtered.loc[all_LAGs_filtered['LAGs'].isna(),'freq.' + sample2].values
NA_group_untreated = all_LAGs_filtered.loc[all_LAGs_filtered['LAGs'].isna(),'freq.' + sample1].values
d_treated = np.log2(all_LAGs_filtered['freq.' + sample2]) - np.log2(gmean(NA_group_treated))
d_untreated = np.log2(all_LAGs_filtered['freq.' + sample1]) - np.log2(gmean(NA_group_untreated))
all_LAGs_filtered['enrichment'] = d_treated - d_untreated

all_gene_LAGs_filtered = add_sense(all_LAGs_filtered.dropna())

all_gene_LAGs_filtered_expanded = expand_multi_LAGs(all_gene_LAGs_filtered)

all_gene_LAGs_filtered_expanded['LAG'] = all_gene_LAGs_filtered_expanded['LAGs']
all_gene_LAGs_filtered_expanded['dir'] = all_gene_LAGs_filtered_expanded['orientation']
all_gene_LAGs_filtered_expanded.to_csv(sample2 + '_vs_' + sample1 + '_expanded_table.txt',sep='\t')

avg_enrichment = pd.DataFrame(all_gene_LAGs_filtered_expanded.groupby(['LAGs','strand'])['enrichment'].mean())
avg_enrichment['n_positions'] = all_gene_LAGs_filtered_expanded.groupby(['LAGs','strand'])['enrichment'].count()
avg_enrichment = avg_enrichment.reset_index()
avg_enrichment = avg_enrichment.loc[avg_enrichment['n_positions']>1]
avg_enrichment = avg_enrichment.sort_values('enrichment')

NA_group = all_LAGs_filtered.loc[all_LAGs_filtered['LAGs'].isna(),'enrichment'].values
## randomly sample non-coding gRNAs to simulate LAGs and get a distribution
print('start sampling')
random.seed(1)
n_position_list = avg_enrichment['n_positions'].unique()
null_list_out = open('source_data/generated/' + sample2 + '_vs_' + sample1 + '_intergenic_null_enrichment.txt','w')
for i in range(0,n_rep):
    for n in n_position_list:
        sampled = random.sample(list(NA_group),n)
        null_list_out.write(str(n) + '\t' + str(np.mean(sampled)) + '\n')
null_list_out.close()

null_list = pd.read_csv('source_data/generated/' + sample2 + '_vs_' + sample1 + '_intergenic_null_enrichment.txt',sep='\t',index_col=0,squeeze=True,header=None)
p_list = []
for i in range(0,len(avg_enrichment)):
    greater = (null_list[avg_enrichment.iloc[i]['n_positions']]>avg_enrichment.iloc[i]['enrichment']).sum()
    less = (null_list[avg_enrichment.iloc[i]['n_positions']]<avg_enrichment.iloc[i]['enrichment']).sum()
    p_list.append(min(greater,less)/(len(null_list[avg_enrichment.iloc[0]['n_positions']])/2))
avg_enrichment['p_val'] = p_list
avg_enrichment.to_csv('source_data/generated/' + sample2 + '_vs_' + sample1 + '_enrichment_table.txt',sep='\t')


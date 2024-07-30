import pandas as pd

#### imports transcript csv file and removed the negative control probes.
#### assigns binary labels to assigned and unassigned cells


def process_data(tf):
    ### filter out negative control probes
    df_neg = tf[tf.feature_name.str.contains('BLANK|Neg', regex=True)].copy()
    df_neg['group'] = 'neg_probes'

    ### filter out transcripts that are genes
    df_genes = tf[~tf['transcript_id'].isin(df_neg.transcript_id)].copy()
    df_genes['group'] = 'gene_probes'
    
    df = pd.concat([df_neg, df_genes], axis=0)

    ### ensure that index for df is equal to original tf
    df.set_index(tf.index, inplace=True)

    ### assign binary labels to assigned and unassigned cells
    df.loc[df.cell_id == 'UNASSIGNED', 'binary'] = 'unassigned'
    df.loc[df.cell_id != 'UNASSIGNED', 'binary'] = 'assigned'
    
    

    return df
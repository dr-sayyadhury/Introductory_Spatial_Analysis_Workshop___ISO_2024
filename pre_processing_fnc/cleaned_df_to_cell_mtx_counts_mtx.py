import pandas as pd

def clean_processed_tf(processed_data, qv=20):
    # Filter for gene probes that are assigned
    gene_assigned = processed_data[(processed_data['binary'] == 'assigned') & (processed_data['group'] == 'gene_probes')].copy()
    # Filter for negative probes that are assigned
    neg_assigned = processed_data[(processed_data['group'] == 'neg_probes') & (processed_data['binary'] == 'assigned')].copy()
    
    # Subset for transcripts with qv > 20
    gene_qv_tf = gene_assigned[gene_assigned['qv'] > qv].copy()

    # Group by cell_id and feature_name, then count the number of transcripts
    gene_qv = gene_qv_tf.groupby(['cell_id', 'feature_name'])['transcript_id'].size().reset_index(name='transcript_count')
    # Pivot table to create a matrix of transcript counts
    gene_mtx = gene_qv.pivot_table(index='cell_id', columns='feature_name', values='transcript_count').fillna(0)
    new_gene_mtx = pd.DataFrame(gene_mtx.values, columns=gene_mtx.columns, index=gene_mtx.index)
    new_gene_mtx.index.name = None
    new_gene_mtx = new_gene_mtx.rename_axis(None, axis=1)

    # Repeat for negative probes
    neg_qv = neg_assigned[neg_assigned['qv'] > qv].copy()
    neg_qv = neg_qv.groupby(['cell_id', 'feature_name'])['transcript_id'].size().reset_index(name='transcript_count')
    neg_mtx = neg_qv.pivot_table(index='cell_id', columns='feature_name', values='transcript_count').fillna(0)
    new_neg_mtx = pd.DataFrame(neg_mtx.values, columns=neg_mtx.columns, index=neg_mtx.index)
    new_neg_mtx.index.name = None    
    new_gene_mtx = new_gene_mtx.rename_axis(None, axis=1)
    
    # Sum the counts across features for each cell
    gene_counts = gene_mtx.sum(axis=1)
    neg_counts = neg_mtx.sum(axis=1)

    df_counts = pd.concat([gene_counts, neg_counts], axis=1)
    df_counts.columns = ['total_counts', 'neg_counts']
    df_counts = df_counts.fillna(0)

    ### calculate centroids
    gene_qv = gene_assigned[gene_assigned['qv'] > qv].copy()
    centroids = gene_qv.groupby('cell_id')[['x_location', 'y_location']].mean().reset_index()
    centroids.columns = ['cell_id', 'centroid_x', 'centroid_y']
    centroids.set_index('cell_id', inplace=True)
    new_centroids = pd.DataFrame(centroids.values, columns=centroids.columns, index=centroids.index)
    new_centroids.index.name = None
    new_centroids = new_centroids.rename_axis(None, axis=1)
    
    return df_counts, gene_qv_tf, new_gene_mtx, new_neg_mtx, centroids

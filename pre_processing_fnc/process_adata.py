import scanpy as sc

def prepare_adata(gene_mtx, df_counts, centroids, min_counts=10, min_cells=3):
    
    adata = sc.AnnData(X=gene_mtx.astype(int).copy())
    df_counts = df_counts[df_counts.index.isin(gene_mtx.index)]

    df_counts = df_counts.reindex(gene_mtx.index)
    adata.obs = df_counts.copy()
    
    centroids = centroids.reindex(adata.obs.index)
    adata.obs[['x_location', 'y_location']] = centroids[['centroid_x', 'centroid_y']].values
    adata.var['genes'] = gene_mtx.columns.astype(str)
    
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)

    return adata


#def normalize_transform_adata(adata, target_sum=1e4, n_comps=21, n_neighbors=10, resolution=0.3):
 
def normalize_transform_adata(adata, args_methods ,args):    
    args_methods = {'SCT': sc.pp.normalize_total, 'log1p': sc.pp.log1p, 'pearson_residuals': sc.experimental.pp.normalize_pearson_residuals}
    args = {'target_sum': 'target_sum', 'n_comps': 'n_comps', 'n_neighbors': 'n_neighbors', 'resolution': 'resolution'}    
    
    for method in args_methods:
        if method in args_methods:
            if method=='SCT':
                args_methods[method](adata)
                adata.layers['SCT'] = adata.X.copy()
                
            elif method=='log1p':
                for size in args['target_sum'] > 0:
                    args_methods[method](adata, target_sum=size)
                    sc.pp.log1p(adata)
                    adata.layers['log1p_'+str(size)] = adata.X.copy()
            
            elif method=='deseq':
                adata.layers['raw_counts'] = adata.X.copy()
                adata.layers['log1p'] = adata.X.copy()
                adata.layers['deseq'] = adata.X.copy()
                
        else:
            print('Method not found')
     

def process_adata(adata, n_comps=21, n_neighbors=10, resolution=0.3):
    sc.pp.pca(adata, n_comps=21)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added='leiden'+str(resolution), resolution=0.3)
    return adata



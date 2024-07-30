import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .qv_labeling import fov_plot  # Import fov_plot

### Example Usage in Another Function:
def plot_gene_and_neg_transcripts(processed_data, type='molecule'):
    gene_assigned = processed_data[(processed_data['binary'] == 'assigned') & (processed_data['group'] == 'gene_probes')].copy()
    neg_assigned = processed_data[(processed_data['group'] == 'neg_probes') & (processed_data['binary'] == 'assigned')].copy()

    if type == 'molecule':
        gene_assigned = gene_assigned.sample(n=neg_assigned.shape[0]*2, random_state=42)
        data_list = [gene_assigned, neg_assigned]

    elif type == 'mean_fov':
        grouped_neg = neg_assigned.groupby('fov_name')
        neg_mean = []
        for fov, group in grouped_neg:
            if group.shape[0] > 0:
                x = group['x_location'].values
                y = group['y_location'].values
                centroid_x = (x.min() + x.max()) / 2
                centroid_y = (y.min() + y.max()) / 2
                qv_avg = group['qv'].mean()
                neg_mean.append({'fov_name': fov, 'x_location': centroid_x, 'y_location': centroid_y, 'qv': qv_avg})
        
        neg_mean = pd.DataFrame(neg_mean)
        
        grouped_gene = gene_assigned.groupby('fov_name')
        gene_mean = []
        for fov, group in grouped_gene:
            if group.shape[0] > 0:
                x = group['x_location'].values
                y = group['y_location'].values
                centroid_x = (x.min() + x.max()) / 2
                centroid_y = (y.min() + y.max()) / 2
                qv_avg = group['qv'].mean()
                gene_mean.append({'fov_name': fov, 'x_location': centroid_x, 'y_location': centroid_y, 'qv': qv_avg})

        gene_mean = pd.DataFrame(gene_mean)
        data_list = [gene_mean, neg_mean]
        
    # Create a figure and axis
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))

    for i, df in enumerate(data_list):
        df = df.sort_values('qv', ascending=True)
        sns.scatterplot(x='x_location', y='y_location', data=df, hue='qv', palette='viridis', s=1, ax=axs[i], legend=False)
        axs[i].set_xlabel('')
        axs[i].set_ylabel('')
        sns.despine()

        norm = plt.Normalize(vmin=0, vmax=40)
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])

        cbar = fig.colorbar(sm, ax=axs[i], ticks=[0, 20, 40])
        cbar.set_label('Quality value', labelpad=15)
        
        fov_plot(processed_data, plot_qv=False, ax=axs[i])
        axs[i].set_title('Gene probes' if i == 0 else 'Negative control probes')
        
    plt.show()


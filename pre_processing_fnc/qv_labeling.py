from matplotlib import colors as mcolors
import matplotlib.pyplot as plt

def fov_plot(processed_data, plot_qv=True, ax=None, identifier=['gene_probes', 'neg_probes']):
    if ax is None:
        fig, ax = plt.subplots()  # Create figure and axis internally if not provided

    if plot_qv:
        qv_min = 0
        qv_max = 40
        cmap = plt.get_cmap('viridis')
        norm = mcolors.Normalize(vmin=qv_min, vmax=qv_max)
    else:
        cmap = None
        norm = None

    grouped = processed_data.groupby('fov_name')
    
    for fov, group in grouped:
        if group.shape[0] > 0:
            x = group['x_location'].values
            y = group['y_location'].values
            xy00 = (x.min(), y.min())
            xy01 = (x.min(), y.max())
            xy10 = (x.max(), y.min())
            xy11 = (x.max(), y.max())
            xy = [xy00, xy01, xy11, xy10, xy00]

            group_assigned = group[group['binary'] == 'assigned']
        
            if plot_qv:
                qv_avg = None
                if (group_assigned['group'] == identifier).any():
                    qv_avg = group_assigned.loc[group_assigned['group'] == identifier, 'qv'].mean()
                elif (group_assigned['group'] == identifier).any():
                    qv_avg = group_assigned.loc[group['group'] == identifier, 'qv'].mean()
                
                if qv_avg is not None:
                    color = cmap(norm(qv_avg))
                    alpha = 0.5
                    ax.fill(*zip(*xy), color=color, alpha=alpha, edgecolor=color, linewidth=0)
                else:
                    color = 'none'
                    alpha = 0
                    ax.fill(*zip(*xy), color=color, alpha=alpha, edgecolor=color, linewidth=0)
            else:
                color = 'none'
                alpha = 0
                ax.fill(*zip(*xy), color=color, alpha=alpha, edgecolor=color, linewidth=0)

            ax.plot(*zip(*xy), color='black', linewidth=0.5)

            centroid_x = (x.min() + x.max()) / 2
            centroid_y = (y.min() + y.max()) / 2

            ax.text(centroid_x, centroid_y, fov, fontsize=8, ha='center', va='center_baseline', color='black')
    
    if plot_qv:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = ax.figure.colorbar(sm, ax=ax)
        cbar.set_label('Quality value (qv)')
    
    if ax is None:
        plt.show()  # Only show plot if created internally

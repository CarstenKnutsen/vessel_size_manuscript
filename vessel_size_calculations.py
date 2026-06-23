
'''
Goal:Create figures for vessel size paper
conda_env:vessel_size
'''
import scanpy as sc
import os
import pandas as pd 
import numpy as np
import seaborn as sns
from functions import compare_obs_values_within_groups_to_excel
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
from scipy.stats import spearmanr
import palantir
from sklearn.preprocessing import MinMaxScaler

adata_name='venous_ec'
figures_base = "data/figures/vessel_size_calculations"
data = "data/single_cell_files/scanpy_files"

os.makedirs(figures_base, exist_ok=True)
sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
sc.settings.figdir = figures_base
sc.settings.autoshow = False

sns.set_style('white', rc={
    'xtick.bottom': True,
    'ytick.left': True,
})
plt.rcParams["font.family"] = "Arial"
size=30



def correlate_genes_with_pseudotime(adata, layer=None, method='spearman',pseudotime='dpt_pseudotime'):
    """
    Correlates all genes with pseudotime in an AnnData object.

    Parameters:
    - adata: AnnData object with pseudotime in `adata.obs['pseudotime']`
    - layer: (Optional) Layer to use instead of adata.X (e.g., 'log1p', 'counts')
    - method: Correlation method, either 'spearman' (default) or 'pearson'

    Returns:
    - pandas DataFrame with genes as index and columns: ['correlation']
    """
    if pseudotime not in adata.obs:
        raise ValueError("Pseudotime must be stored in adata.obs['pseudotime'].")

    # Get expression matrix
    X = adata.X if layer is None else adata.layers[layer]
    if not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(X.toarray() if hasattr(X, "toarray") else X,
                         index=adata.obs_names, columns=adata.var_names)

    # Extract pseudotime
    pseudotime = adata.obs[pseudotime]

    # Run correlation
    results = []
    for gene in X.columns:
        if method in ['spearman','pearson']:
            corr = X[gene].corr(pseudotime,method=method)  # Pearson p-value not computed here
        else:
            raise ValueError("Method must be 'spearman' or 'pearson'.")
        results.append((gene, corr))

    result_df = pd.DataFrame(results, columns=['gene', 'correlation']).set_index('gene')
    return result_df.sort_values('correlation', ascending=False)

scaler = MinMaxScaler()
def normalize_dataframe(df, feature_range = (0, 1)):
    # Initialize the MinMaxScaler
    scaler = MinMaxScaler(feature_range=feature_range)
    # Fit the scaler on the data and transform each column
    df_normalized = pd.DataFrame(scaler.fit_transform(df), index=df.index, columns=df.columns)
    return df_normalized


if __name__ == '__main__':
    adata_fns = {'all':'venous_ec_vessel_size_velocity.gz.h5ad',
                 'Normoxia':'venous_ec_vessel_size_Normoxia_velocity.gz.h5ad',
                 'Hyperoxia':'venous_ec_vessel_size_Hyperoxia_velocity.gz.h5ad'
                 }
    for run in adata_fns:
        fn = adata_fns[run]
        adata = sc.read(f"{data}/{fn}")
        if run =='all':
            adata.obs['Treatment'] = pd.Categorical(adata.obs['Treatment'], categories=['Normoxia', 'Hyperoxia'])
        del_layers = ['MAGIC_imputed_data', 'Ms', 'Mu', 'ambiguous', 'cp10k', 'fit_t', 'fit_tau', 'fit_tau_', 'log1p_cc_regress','matrix', 'raw', 'spliced', 'unspliced', 'variance_velocity', 'velocity', 'velocity_u']
        for x in del_layers:
            del adata.layers[x]
        if run != 'all':
            figures = f"{figures_base}/{run}"
            os.makedirs(figures, exist_ok=True)
            sc.settings.figdir = figures
        else:
            figures=figures_base


        ## correlate each ct with pseudotime
        corr_dfs = {}
        for ct in ['Arterial EC','Venous EC']:
            ct_adata = adata[adata.obs['Cell Subtype_no_cc']==ct]
            df = correlate_genes_with_pseudotime(ct_adata,layer='log1p',method='pearson',pseudotime='palantir_pseudotime')
            df = df.dropna(how='all')
            df.to_csv(f'{figures}/pseudotime_correlation_{ct}.csv')
            corr_dfs[ct]=df

        # find common correlations and find overlap
        if run == 'all':
            for top_n_genes in [50,150]:
                arterial_large_genes = corr_dfs['Arterial EC'].head(top_n_genes).index.tolist()
                venous_large_genes = corr_dfs['Venous EC'].head(top_n_genes).index.tolist()
                arterial_small_genes = corr_dfs['Arterial EC'].tail(top_n_genes).index.tolist()[::-1]
                venous_small_genes = corr_dfs['Venous EC'].tail(top_n_genes).index.tolist()[::-1]

                ##Create the Venn diagrams
                ## large top 50
                venn = venn2([set(arterial_large_genes), set(venous_large_genes)],
                             set_labels=('Arterial', 'Venous'),
                             set_colors=('#4A90E2', '#E35D6A'),
                             alpha=0.7)
                for text in venn.set_labels:
                    text.set_fontsize(12)
                for text in venn.subset_labels:
                    if text:
                        text.set_fontsize(12)
                plt.title(f"Top {top_n_genes} genes positively correlated with pseudotime")
                plt.savefig(f'{figures}/venn_diagram_large_{top_n_genes}.png',dpi=300,bbox_inches='tight')
                plt.close()

                ## Create the Venn diagrams
                ## small top 50
                venn = venn2([set(arterial_small_genes), set(venous_small_genes)],
                             set_labels=('Arterial', 'Venous'),
                             set_colors=('#4A90E2', '#E35D6A'),
                             alpha=0.7)
                for text in venn.set_labels:
                    text.set_fontsize(12)
                for text in venn.subset_labels:
                    if text:
                        text.set_fontsize(12)
                plt.title(f"Top {top_n_genes} genes negatively correlated with pseudotime")
                plt.savefig(f'{figures}/venn_diagram_small_{top_n_genes}.png',dpi=300,bbox_inches='tight')
                plt.close()

                ## calculate  vessel size
                large_genes = [x for x in arterial_large_genes if x in venous_large_genes]
                small_genes = [x for x in arterial_small_genes if x in venous_small_genes]
                if top_n_genes==50:
                    adata.uns['large_genes'] = large_genes
                    adata.uns['small_genes'] = small_genes
                else:
                    adata.uns[f'large_genes_{top_n_genes}'] = large_genes
                    adata.uns[f'small_genes_{top_n_genes}'] = small_genes

                genes = large_genes + small_genes
                sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))

                sc.pl.umap(adata,color=genes,hspace=0.5,save=f'size_genes_{top_n_genes}.png')
                corr_df_50 = pd.DataFrame([small_genes,large_genes],index=None).T
                corr_df_50.columns = ['small_genes','large_genes']
                corr_df_50.to_csv(f'{figures}/top{top_n_genes}_shared_corr_df.csv')
                sc.tl.score_genes(adata,large_genes,score_name='large_score')
                sc.tl.score_genes(adata,small_genes,score_name='small_score')
                if top_n_genes==50:
                    label = 'Vessel size'
                else:
                    label = f'Vessel size {top_n_genes}'
                adata.obs[f'{label} score'] = adata.obs['large_score'] - adata.obs['small_score']
                adata.obs[f'{label} score'] = normalize_dataframe(adata.obs[[f'{label} score']])
                adata.obs[f'{label} category'] = pd.cut(adata.obs[f'{label} score'], bins=4,labels=['capillary','small','medium','large'])
                genes = large_genes + small_genes
                fig = palantir.plot.plot_gene_trend_heatmaps(adata, genes, cmap='viridis')
                fig.tight_layout()
                fig.savefig(f'{figures}/palantir_heatmap_gene_trends_{top_n_genes}.png')
                plt.close()
        print(adata)
        sc.pl.umap(adata,
                                  color='Cell Subtype_no_cc',
                                  na_in_legend=False,
                                  size=size,
                                  frameon=False,
                                  title="",
                                  # legend_loc=None,
                                  legend_fontsize=10,
                                  legend_fontoutline=0,
                                  save='vascular_celltypes.png'
                                  )
        sc.pl.umap(adata,
                   color='Treatment',
                   na_in_legend=False,
                   size=size,
                   frameon=False,
                   title="",
                   # legend_loc=None,
                   legend_fontsize=10,
                   legend_fontoutline=0,
                   save='vascular_treatment.png'
                   )
        sc.pl.umap(adata,color=['Vessel size score'],cmap='Oranges',frameon=False,size=size,save='vessel_size_score.png')
        sc.pl.umap(adata,color=['palantir_pseudotime'],cmap='viridis',frameon=False,size=size,save='palantir_pseudotime.png')

        # set ct_s_pos to delineate locations
        adata_cap_ident = adata[adata.obs['Cell Subtype_no_cc'] == 'Cap1']
        adata_cap_ident.obs['Capillary side'] = ['Cap1-Vein' if x == True else 'Cap1-Artery' for x in
                                                 adata_cap_ident.obsm['branch_masks']['Venous EC']]
        adata.obs['Capillary side'] = adata_cap_ident.obs['Capillary side']
        ls = []
        for x,y,z in zip(adata.obs['Cell Subtype_no_cc'],adata.obs['Vessel size category'],adata.obs['Capillary side']):
            if x =='Cap1':
                ls.append(z)
            else:
                if y =='capillary':
                    if x == 'Arterial EC':
                        ls.append(f'Cap1-Artery')
                        continue
                    else:
                        ls.append(f'Cap1-Vein')
                        continue
                if x =='Arterial EC':
                    x = 'PAEC'
                else:
                    x ='PVEC'
                ls.append(f'{x} {y[0].upper()}')
        adata.obs['ct_s_pos'] = ls

        # generate some marker gene lists
        sc.tl.rank_genes_groups(adata, groupby="ct_s_pos", method='wilcoxon', pts=True)
        with pd.ExcelWriter(f"{figures}/ct_size_pos_markers.xlsx", engine="xlsxwriter") as writer:
            for ct in sorted(adata.obs["ct_s_pos"].unique()):
                df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups", group=ct)
                df.set_index("names")
                df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
                df.sort_values('pct_difference',ascending=False).to_excel(writer, sheet_name=f"{ct} v rest"[:31])


        sc.tl.rank_genes_groups(adata, groupby="Vessel size category", method='wilcoxon', pts=True)
        with pd.ExcelWriter(f"{figures}/size_markers.xlsx", engine="xlsxwriter") as writer:
            for ct in sorted(adata.obs["Vessel size category"].unique()):
                df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups", group=ct)
                df.set_index("names")
                df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
                df.sort_values('pct_difference',ascending=False).to_excel(writer, sheet_name=f"{ct} v rest"[:31])
        df = adata.obs[['palantir_pseudotime','Vessel size score','Vessel size 150 score','ct_s_pos']]
        df.to_csv(f'{figures}/select_obs.csv')
        if run =='all':
            adata.write(f'{data}/{adata_name}_vessel_size_plot.gz.h5ad',compression='gzip')
            fig, axs = plt.subplots(1, 2, figsize=(4, 2))
            axs = axs.ravel()
            plt.subplots_adjust(wspace=0.75)
            delta = df['Vessel size score'] - df['Vessel size 150 score']
            sns.regplot(x=df['Vessel size score'], y=df['Vessel size 150 score'], scatter_kws={'s': 3, 'linewidths': 0},
                        line_kws={'color': 'red'}, ci=None, ax=axs[0])
            axs[0].set_xlabel('Vessel size score\50')
            axs[0].set_ylabel('Vessel size score\n150')
            r, p = sp.stats.pearsonr(df['Vessel size score'], df['Vessel size 150 score'])
            axs[0].text(0.5, 1.05, 'r = {:.2f}'.format(r, p), fontsize=10,
                        transform=axs[0].transAxes)
            sns.violinplot(delta, ax=axs[1])
            axs[1].set_ylabel('Vessel size score delta')
            fig.savefig(f'{figures}/vessel_size_calculation_top_n_genes_comparison.png', bbox_inches='tight', dpi=300)
        else:
            adata.write(f'{data}/{adata_name}_vessel_size_{run}_plot.gz.h5ad', compression='gzip')

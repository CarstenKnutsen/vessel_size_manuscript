
'''
Goal:Integrate all vessel size datasets together
'''
import scanpy as sc
import scanpy.external as sce
import os 
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import palantir
import gprofiler
import anndata
from gprofiler import GProfiler
from upsetplot import from_contents, UpSet
import scipy as sp

figures = "data/figures/figures/"
figures_out = f'{figures}/integrate_all'
os.makedirs(figures_out, exist_ok=True)
sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
sc.settings.figdir = figures_out
sns.set_style('white', rc={
    'xtick.bottom': True,
    'ytick.left': True,
})
plt.rcParams["font.family"] = "Arial"
size=5
sc.settings.autoshow = False
def convert_mouse_to_human(adata):
    """Convert var_names in an AnnData object from mouse to human gene names using gProfiler."""
    gp = GProfiler(return_dataframe=True)

    # Convert mouse genes to human
    conversion_df = gp.orth(organism='mmusculus',query=adata.var_names.tolist(), target='hsapiens')

    # Extract mapping
    mouse_to_human = dict(zip(conversion_df['incoming'], conversion_df['name']))
    # Map `var_names`
    adata.var_names = adata.var_names.map(lambda x: mouse_to_human.get(x) if mouse_to_human.get(x) not in [None, 'N/A'] else x)
    adata.var_names = adata.var_names.str.upper() #Capitalize any unmapped genes
    # Remove duplicates (if any)
    adata.var_names_make_unique()
    return adata
if __name__ == '__main__':
    adata_dict = {}
    adata = sc.read(f'{figures}/hurskainen_2021/vessel_size.gz.h5ad')
    adata.x = adata.layers['log1p'].copy()
    del adata.layers
    adata.var['mouse_genes'] = adata.var_names.copy()
    adata = convert_mouse_to_human(adata)
    adata.obs['Species'] = 'Mouse'
    adata.obs['Age'] = 'Neonatal'
    adata_dict['Hurskainen 2021'] = adata

    adata = sc.read(f'{figures}/tabula_muris_senis/vessel_size.gz.h5ad')
    del adata.layers
    adata.var['mouse_genes'] = adata.var_names.copy()
    adata = convert_mouse_to_human(adata)
    del adata.raw
    adata.obs['Species'] = 'Mouse'
    adata.obs['Age'] = 'Adult'
    adata_dict['Tabula Muris Senis 2020'] = adata

    adata = sc.read(f"data/single_cell_files/scanpy_files/venous_ec_vessel_size_plot.gz.h5ad")
    adata.X = adata.layers['log1p']
    del adata.layers
    adata.var['mouse_genes'] = adata.var_names.copy()
    adata = convert_mouse_to_human(adata)
    adata.obs['Cell Subtype'] = adata.obs['Cell Subtype_no_cc'].copy()
    adata.obs['Species'] = 'Mouse'
    adata.obs['Age'] = 'Neonatal'
    adata_dict['Sveiven 2026'] = adata
    adata = sc.read(f'{figures}/bhattacharya_2024/vessel_size.gz.h5ad')
    del adata.layers
    adata.obs['Species'] = 'Human'
    adata.obs['Age'] = 'Neonatal'
    adata_dict['Bhattacharya 2024'] = adata
    adata = sc.read(f'{figures}/tabula_sapiens/vessel_size.gz.h5ad')
    del adata.layers
    adata.obs['Species'] = 'Human'
    adata.obs['Age'] = 'Adult'
    adata_dict['Tabula Sapiens v2 2025'] = adata
    adata = sc.read(f'{figures}/lungmap_bpd_2024/vessel_size.gz.h5ad')
    del adata.layers
    adata_dict['LungMAP 2024'] = adata
    adata.obs['Species'] = 'Human'
    adata.obs['Age'] = 'Neonatal'
    adata = sc.read(f'{figures}/zhao_2024/vessel_size.gz.h5ad')
    adata.x = adata.layers['log1p'].copy()
    del adata.layers
    adata.var['mouse_genes'] = adata.var_names.copy()
    adata = convert_mouse_to_human(adata)
    adata.obs['Species'] = 'Mouse'
    adata.obs['Age'] = 'Adult'
    adata_dict['Zhao 2024'] = adata
    ## Create upset plot to show vessel size signature overlap
    gene_dict = {}
    gp = GProfiler(return_dataframe=True)
    for ds in adata_dict:
        adata = adata_dict[ds]
        print(ds)
        gene_ls = adata.uns['large_genes'].tolist() + adata.uns['small_genes'].tolist()
        if 'Mouse' in adata.obs['Species'].unique():
            conversion_df = gp.orth(organism='mmusculus', query=gene_ls, target='hsapiens')
            conversion_df.loc[conversion_df["name"] == "N/A", "name"] = conversion_df.loc[
                conversion_df["name"] == "N/A", "incoming"].str.upper()
            gene_ls = conversion_df.name.values.tolist()
        gene_dict[ds] = gene_ls
    upset = from_contents(gene_dict)
    membership_df = upset.reset_index()
    membership_df = membership_df.set_index('id')
    membership_df['# of datasets used for scoring'] = membership_df.sum(axis=1)
    membership_df.sort_values('# of datasets used for scoring', ascending=False).to_csv(f"{figures_out}/membership_df.csv")
    print(upset)
    ax_dict = UpSet(upset, sort_by='cardinality',subset_size="count").plot()
    plt.savefig(f'{figures_out}/upset_plot_shared_genes_signature.png', dpi=300, bbox_inches='tight')
    plt.close()
    ## make overlap barplot
    membership_df = upset.reset_index()
    membership_df = membership_df.set_index('id')
    overlap_counts = {}
    total_true = membership_df.sum(axis=0)

    for col in membership_df.columns:
        # True in any OTHER dataset
        other_true = membership_df.drop(columns=col).any(axis=1)

        # True in this dataset AND another dataset
        overlap = membership_df[col] & other_true
        unique = membership_df[col] & ~other_true

        overlap_counts[col] = overlap.sum()

    overlap_counts = pd.Series(overlap_counts)

    plot_df = pd.DataFrame({
        "Dataset": membership_df.columns,
        "Total genes used": total_true.values,
        "Overlap with 1+ dataset": overlap_counts.values
    })

    plot_long = plot_df.melt(
        id_vars="Dataset",
        var_name="Category",
        value_name="Count"
    )

    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    plot_long = plot_long.sort_values('Dataset')
    ax = sns.barplot(
        data=plot_long,
        x="Dataset",
        y="Count",
        hue="Category",
        dodge=False,
    )
    print(plot_df)
    print(plot_long)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    sns.despine()
    plt.legend(frameon=False)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(0.55, 1))
    ax.set_ylabel('# genes used\nin vessel size score')
    ax.set_xlabel('')
    fig.savefig(f'{figures_out}/barplot_overlap_dataset.png', dpi=300, bbox_inches='tight')

    ## Make same plots of genes to compare across datasets
    gene_ls_dict = {'shared_size': [
        'Ccdc85a',
        'Glp1r',
        'Kit',
        'Sox4', 'Nrp1', 'Ifitm3',
        'Ptprr', 'Adgrg6', 'Foxo1',
        'Mgp', 'Eln', 'Nr4a2',
    ][::-1],
                    'size_axis': ['Dkk2', 'Sox6',
                                  'Lama3', 'Bdkrb2',
                                  'Stc1',
                                  'Adam23', 'Ntrk2',
                                  'Emid1', 'Chrm2',
                                  'Chrm3', 'Rarb',
                                  'Gria3',
                                  'Ptger3', 'Moxd1',
                                  ],
                    'umap_genes':['Eln','Mgp','Mecom','Col4a1','Col4a2']
                    }
    for ds in adata_dict:
        if ds == 'Sveiven 2026':
            continue
        print(ds)
        adata = adata_dict[ds]
        for x in gene_ls_dict:
            gene_ls = gene_ls_dict[x].copy()
            if 'Human' in adata.obs['Species'].unique():
                gp = GProfiler(return_dataframe=True)
                conversion_df = gp.orth(organism='mmusculus', query=gene_ls, target='hsapiens')
                conversion_df.loc[conversion_df["name"] == "N/A", "name"] = conversion_df.loc[
                    conversion_df["name"] == "N/A", "incoming"].str.upper()
                gene_ls = conversion_df.name.values.tolist()
                gene_ls = [x for x in gene_ls if x in adata.var_names]
                gene_symbols = None
            else:
                gene_ls = [x for x in gene_ls if x in adata.var['mouse_genes'].values]
                gene_symbols = 'mouse_genes'

            sc.pl.DotPlot(adata, gene_ls, groupby='ct_s',
                          gene_symbols=gene_symbols,
                          title=ds,
                          categories_order=['PAEC L', 'PAEC M', 'PAEC S', 'Cap1', 'PVEC S', 'PVEC M', 'PVEC L'],
                          standard_scale='var').style(cmap='viridis').savefig(
                f'{figures_out}/{x}_{ds}_dotplot_size_axis.png', dpi=300,
                bbox_inches='tight')
            if x =='umap_genes':
                for gene in gene_ls:
                    sc.pl.umap(adata, color=gene, gene_symbols=gene_symbols,cmap='viridis', frameon=False, size=size, save=f'_{gene}_{ds}.png')

    adata = anndata.concat(adata_dict.values(), join='inner', label='Dataset', keys=adata_dict.keys(), index_unique=None)
    del adata.obsm
    adata.obs_names_make_unique()
    nindex = np.random.permutation(adata.obs.index)
    adata=adata[nindex,:]
    ## read in cell cycle genes
    cell_cycle_genes = [x.strip() for x in open('/home/carsten/alvira_bioinformatics/venous_ec_scrnaseq/data/outside_data/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    adata.layers['log1p'] = adata.X.copy()
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    adata.uns['Cell Subtype_colors']= ['#4A90E2','#9B59B6','#E35D6A']
    sc.pp.highly_variable_genes(adata, batch_key='Dataset',n_top_genes=200)
    sc.pp.pca(adata,use_highly_variable=True)
    sce.pp.bbknn(adata,batch_key='Dataset',neighbors_within_batch=2,set_op_mix_ratio=5)
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata,min_dist=0.5)
    adata.X = adata.layers['log1p'].copy()
    sc.tl.rank_genes_groups(adata,'leiden',method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(adata,dendrogram=False)
    sc.pl.umap(adata,color=['leiden'],cmap='Oranges',size=size,frameon=False,save='_all_leiden.png')
    sc.pl.umap(adata,color=['Vessel size score'],cmap='Oranges',size=size,frameon=False,save='_all_vessel_size_score.png')
    sc.pl.umap(adata,color=['Vessel size category'],cmap='viridis',size=size,frameon=False,save='_allvessel_size_category.png')
    sc.pl.umap(adata,color=['Cell Subtype'],cmap='viridis',size=size,legend_fontsize=10, legend_fontoutline=1,frameon=False,save='_all_cellsubtype.png')
    sc.pl.umap(adata,color=['MGP'],cmap='viridis',frameon=False,size=size,save='_all_mgp.png')
    sc.pl.umap(adata,color=['COL4A1'],cmap='viridis',frameon=False,size=size,save='_all_col4a1.png')
    sc.pl.umap(adata,color=['COL4A2'],cmap='viridis',frameon=False,size=size,save='_all_col4a2.png')
    sc.pl.umap(adata,color=['ELN'],cmap='viridis',frameon=False,size=size,save='_all_eln.png')
    sc.pl.umap(adata,color=['FBLN2'],cmap='viridis',frameon=False,size=size,save='_all_fbln2.png')
    sc.pl.umap(adata,color=['HEY1'],cmap='viridis',frameon=False,size=size,save='_all_hey1.png')
    sc.pl.umap(adata,color=['NR2F2'],cmap='viridis',frameon=False,size=size,save='_all_nr2f2.png')
    sc.pl.umap(adata,color=['Dataset'],cmap='viridis',frameon=False,size=size,save='_all_dataset.png')
    sc.pl.umap(adata,color=['MKI67'],cmap='viridis',frameon=False,size=size,save='_all_mki67.png')
    sc.pl.umap(adata,color=['Age'],palette=['purple','yellow'],frameon=False,size=size,save='_all_age.png')
    sc.pl.umap(adata,color=['Species'],cmap='viridis',frameon=False,size=size,save='_all_species.png')
    sc.pl.umap(adata,color=['Dataset'],cmap='viridis',frameon=False,size=size,save='_all_dataset.png')
    sc.pl.umap(adata,color='Vessel size category')
    for dataset in adata.obs['Dataset'].unique():
        sc.pl.umap(adata.copy(),color=['Dataset'],cmap='viridis',mask_obs=adata.obs['Dataset']==dataset,title='',legend_loc=None,frameon=False,size=size,save=F'_{dataset}_dataset.png')
        sc.pl.umap(adata.copy(),color='Vessel size category',mask_obs=adata.obs['Dataset']==dataset,title='',legend_loc=None,frameon=False,size=size,save=F'_{dataset}_dataset_vessel_size.png')
        sc.pl.umap(adata.copy(),color='Vessel size score',mask_obs=adata.obs['Dataset']==dataset,cmap='Oranges',title='',colorbar_loc=None,frameon=False,size=size,save=F'_{dataset}_dataset_vessel_size_score.png')
    ## create scatter and violin to show vessel size score correlation between neighbors of datasets
    obs_key = "Vessel size score"
    dataset_key = "Dataset"
    vals = adata.obs[obs_key].values
    datasets = adata.obs[dataset_key].values
    conn = adata.obsp["connectivities"]
    neighbor_median = np.zeros(adata.n_obs)
    for i in range(adata.n_obs):
        neighbors = conn[i].indices  # neighbor indices
        valid_neighbors = neighbors[datasets[neighbors] != datasets[i]]
        if len(valid_neighbors) > 0:
            neighbor_median[i] = np.median(vals[valid_neighbors])
        else:
            neighbor_median[i] = np.nan  # no cross-dataset neighbors
    adata.obs[f"{obs_key}_neighbor_median_cross_dataset"] = neighbor_median
    adata.write(f'{figures_out}/vessel_size_all_datasets.gz.h5ad',compression='gzip')
    ## make plots showing correlation
    df = adata.obs[['Vessel size score', 'Vessel size score_neighbor_median_cross_dataset', 'Dataset']]
    df = df.dropna()
    df['bins'] = pd.cut(df['Vessel size score'], bins=10)
    joint_plot = sns.jointplot(x=df['Vessel size score'], y=df['Vessel size score_neighbor_median_cross_dataset'],
                               kind='reg',
                               scatter_kws={'s': 0.1, 'linewidths': 0},
                               line_kws={'color': 'red'}, height=2)
    joint_plot.set_axis_labels("Vessel size score", "Vessel size score\nneighbor median")
    r, p = sp.stats.pearsonr(df['Vessel size score'], df['Vessel size score_neighbor_median_cross_dataset'])
    joint_plot.ax_marg_x.text(0.5, 1.05, 'r = {:.2f}'.format(r, p), fontsize=10,
                              transform=joint_plot.ax_marg_x.transAxes)
    joint_plot.ax_joint.set_ylim([0, 1])
    joint_plot.ax_joint.set_yticks([0, 0.5, 1])
    joint_plot.savefig(f'{figures_out}/jointplot_neighbors.png', dpi=300, bbox_inches='tight')
    plt.close
    fig, ax = plt.subplots(1, 1, figsize=(2, 2))

    sns.violinplot(df, x='bins', y='Vessel size score_neighbor_median_cross_dataset', ax=ax)
    ax.set_xticklabels([f'{str(x)[:3]}-{str(x + 0.1)[:3]}' for x in np.arange(0, 1, 0.1)]
                       , rotation=90)
    ax.set_xlabel('Vessel size score')
    ax.set_ylabel('Vessel size score\nneighbor median')
    fig.savefig(f'{figures_out}/violinplot_neighbors.png', dpi=300, bbox_inches='tight')

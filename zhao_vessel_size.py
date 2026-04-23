
'''
Goal:Check tabula sapiens lung for vessel_size_gradient
conda_env:vessel_size

'''

import scanpy as sc
import scanpy.external as sce
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import palantir
import cellrank as cr
import scvelo as scv
from functions import palantir_psuedotime_routine_external_datasets,find_gene_overlap_in_pseudotimes
from scipy.stats import median_abs_deviation

adata_name='zhao_2024'
figures = "data/figures/figures/zhao_2024"
data = "data/single_cell_files/scanpy_files"

os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
sc.settings.figdir = figures
sns.set_style('white', rc={
    'xtick.bottom': True,
    'ytick.left': True,
})
plt.rcParams["font.family"] = "Arial"
size=15
sc.settings.autoshow = False


if __name__ == '__main__':
    adata = sc.read_10x_h5(
        '/home/carsten/alvira_bioinformatics/data/external_datasets/zhao2024/GSE201631_RAW/GSM6068821_D0_filtered_feature_bc_matrix.h5')
    adata.var_names_make_unique()
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["mt"] = adata.var_names.str.startswith(("mt"))
    adata.var["hb"] = adata.var_names.str.contains(("^Hb[^(p)]"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        expr_type="umis",
        percent_top=[20],
        log1p=True,
        inplace=True,
    )


    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
                np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier


    adata.obs["outlier"] = (is_outlier(adata, "log1p_total_umis", 5) | is_outlier(adata, "log1p_n_genes_by_umis",
                                                                                  5)) | is_outlier(adata,
                                                                                                   "pct_umis_in_top_20_genes",
                                                                                                   5)
    adata.obs["mt_outlier"] = adata.obs["pct_umis_mt"] > 30

    print(adata.obs.mt_outlier.value_counts())
    print(adata.obs.outlier.value_counts())
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(adata.obs.total_umis.median())
    print(adata.obs.n_genes_by_umis.median())
    sc.pp.normalize_total(adata, key_added=None, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata, mask_var="highly_variable")
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, key_added="leiden", resolution=0.5)
    sc.tl.umap(adata, min_dist=0.5)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", pts=True)
    sc.pl.rank_genes_groups_dotplot(adata, save='initial_leiden_markers.png')
    with pd.ExcelWriter(
            f"{figures}/initial_leiden_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs["leiden"].cat.categories:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups", group=ct)
            df.set_index("names")
            df["pct_difference"] = df["pct_nz_group"] - df["pct_nz_reference"]
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
    sc.pl.umap(adata, color=[
                             'leiden','Cdh5',
                             'Ptprc',
                             'Epcam',
                             'Col1a1',
                             'Col3a1',
                             'Gja5',
                             'Kit',
                             'Slc6a2',
                             'Tbx2',
                             'Car4',
                             'Col15a1',
                             'Fabp4',
                             'Prox1',
                             'Mki67',
                             'pct_umis_ribo', 'log1p_n_genes_by_umis'], hspace=0.5,wspace=0.5,save='_zhao2024_initial_genes.png')
    leiden_dict = {
        "0": "Cap1",
        "1": "Cap1",
        "2": "Cap2",
        "3": "Venous EC",
        "4": "Arterial EC",
        "5": "Cap1",
        "6": "Lymphatic EC",
        "7": "Immune",
        "8": "Immune",
        "9": "Cap2",
        "10": "Mesenchymal doublets",
        "11": "low-quality low umis genes",
        "12": "Cap1",
    }
    adata.obs['Cell Subtype'] = adata.obs['leiden'].map(leiden_dict)
    vessel_cts = ['Arterial EC', 'Cap1', 'Venous EC']
    adata_vessel = adata[adata.obs['Cell Subtype'].isin(vessel_cts)]
    sc.pp.highly_variable_genes(adata_vessel, n_top_genes=500)
    sc.pp.pca(adata_vessel, mask_var="highly_variable")
    sc.pp.neighbors(adata_vessel)
    sc.tl.leiden(adata_vessel, key_added='leiden_vessel', resolution=0.5)
    sc.tl.umap(adata_vessel, min_dist=0.5)
    sc.pl.umap(adata_vessel, color=['Cdh5',
                                    'leiden',
                                    'leiden_vessel',
                                    "Cell Subtype",
                                    'Eln',
                                    'Col4a1',
                                    'Ptprc',
                                    'Epcam',
                                    'Col1a1',
                                    'Col3a1',
                                    'Gja5',
                                    'Kit',
                                    'Slc6a2',
                                    'Tbx2',
                                    'Car4',
                                    'Col15a1',
                                    'Fabp4',
                                    'Prox1',
                                    'Mki67',
                                    'pct_umis_ribo', 'log1p_n_genes_by_umis'], hspace=0.5,save ='meta_umaps.png')
    plt.close()
    adata_vessel.obs['Cell Subtype'] = adata_vessel.obs['leiden_vessel'].map({'0':'Cap1','1':'Cap1','2':'Venous EC','3':'Arterial EC'})
    adata_vessel.uns['Cell Subtype_colors'] = ['#4A90E2', '#9B59B6', '#E35D6A']
    adata_name = 'zhao2024'
    start_and_end_states = {'root_ct': ('Cap1', 1, 'min'),
                            'terminal_cts': {'Arterial EC': (1, 'max'),
                                             'Venous EC': (1, 'max')
                                             }
                            }
    adata_vessel.obsm['X_pca_harmony'] = adata_vessel.obsm['X_pca'].copy()  # only one library this dataset
    palantir_psuedotime_routine_external_datasets(adata_vessel, start_and_end_states, save_prefix=figures,
                                                  adata_name=adata_name)
    find_gene_overlap_in_pseudotimes(adata_vessel, save_prefix=figures, adata_name=adata_name)
    ls = []
    for x, y in zip(adata_vessel.obs['Cell Subtype'], adata_vessel.obs['Vessel size category']):
        if x == 'Cap1':
            ls.append('Cap1')
        else:
            if y == 'capillary':
                ls.append(f'Cap1')
                continue
            if x == 'Arterial EC':
                x = 'PAEC'
            else:
                x = 'PVEC'
            ls.append(f'{x} {y[0].upper()}')
    adata_vessel.obs['ct_s'] = ls
    gene_ls = [
        'Ccdc85a',
        'Glp1r',
        'Kit',
        'Sox4', 'Nrp1', 'Ifitm3',
        'Ptprr', 'Adgrg6', 'Foxo1',
        'Mgp', 'Eln', 'Nr4a2',
    ][::-1]
    sc.pl.umap(adata_vessel, color=gene_ls, hspace=0.5, wspace=0.3, cmap='viridis', save='size_axis.png')
    sc.pl.MatrixPlot(adata_vessel, gene_ls, groupby='ct_s',
                     categories_order=['PAEC L', 'PAEC M', 'PAEC S', 'Cap1', 'PVEC S', 'PVEC M', 'PVEC L'],
                     standard_scale='var').style(cmap='viridis').swap_axes().savefig(
        f'{figures}/matrixplot_size_axis.png', dpi=300, bbox_inches='tight')
    sc.pl.DotPlot(adata_vessel, gene_ls, groupby='ct_s',
                  categories_order=['PAEC L', 'PAEC M', 'PAEC S', 'Cap1', 'PVEC S', 'PVEC M', 'PVEC L'],
                  standard_scale='var').style(cmap='viridis').savefig(f'{figures}/dotplot_size_axis.png', dpi=300,
                                                                      bbox_inches='tight')
    # In[31]:

    gene_ls = ['Dkk2', 'Sox6',
               'Lama3', 'Bdkrb2',
               # 'Depp1', not present
               'Stc1',
               'Adam23', 'Ntrk2',
               'Emid1', 'Chrm2',
               'Chrm3', 'Rarb',
               # 'Fads2b', not present
               'Gria3',
               'Ptger3', 'Moxd1',
               ]
    sc.pl.MatrixPlot(adata_vessel, gene_ls, groupby='ct_s',
                     categories_order=['PAEC L', 'PAEC M', 'PAEC S', 'Cap1', 'PVEC S', 'PVEC M', 'PVEC L'],
                     standard_scale='var').style(cmap='viridis').swap_axes().savefig(
        f'{figures}/matrixplot_size_markers.png', dpi=300, bbox_inches='tight')
    sc.pl.DotPlot(adata_vessel, gene_ls, groupby='ct_s',
                  categories_order=['PAEC L', 'PAEC M', 'PAEC S', 'Cap1', 'PVEC S', 'PVEC M', 'PVEC L'],
                  standard_scale='var').style(cmap='viridis').savefig(f'{figures}/dotplot_size_markers.png', dpi=300,
                                                                      bbox_inches='tight')




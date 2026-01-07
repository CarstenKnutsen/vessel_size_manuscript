


import scanpy as sc
import os
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import liana as li

adata_name='venous_ec'
figures = "data/figures/figures"
data = "data/single_cell_files/scanpy_files"

os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi_save=300, fontsize=8, figsize=(2,2))
sc.settings.figdir = figures
sns.set_style('white', rc={
    'xtick.bottom': True,
    'ytick.left': True,
})
plt.rcParams["font.family"] = "Arial"
size=15


if __name__ == '__main__':
    adata = sc.read(f"{data}/{adata_name}_vessel_size_plot.gz.h5ad")
    ct_s_pos = adata.obs['ct_s_pos']
    del adata
    adata_all = sc.read(f'{data}/{adata_name}_celltyped_no_cc.gz.h5ad')
    adata_all.obs['Cell Subtype_w_size'] = adata_all.obs['Cell Subtype_no_cc'].astype(str).copy()
    adata_all.obs.loc[ct_s_pos.index, 'Cell Subtype_w_size'] = ct_s_pos.astype(str)
    li.mt.rank_aggregate(adata_all,
                         groupby='Cell Subtype_w_size',
                         resource_name='mouseconsensus',
                         expr_prop=0.1,
                         use_raw=False,
                         verbose=True)
    df = adata_all.uns['liana_res'].copy()
    df.to_csv(f'{figures}/liana_vessel_size_ccc.csv')


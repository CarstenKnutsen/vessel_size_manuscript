
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
from functions import palantir_psuedotime_routine_external_datasets,find_gene_overlap_in_pseudotimes,

adata_name='tabula_muris_senis'
figures = "data/figures/figures/tabula_muris_senis"
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

    adata_tms = sc.read('/home/carsten/alvira_bioinformatics/data/external_datasets/tabula_muris_senis/raw_data/tabula_muris_senis_lung_droplet.h5ad')
    del adata_tms.raw
    sc.pp.normalize_total(adata_tms,target_sum=1e4)
    sc.pp.log1p(adata_tms)
    adata_tms.var_names = adata_tms.var['feature_name'].str.split('_').str[0]
    adata_tms.var_names_make_unique()
    adata_tms.var.index.name = 'gene'
    vessel_cts = ['Artery','Capillary','Vein']
    adata_tms_vessel = adata_tms[adata_tms.obs['free_annotation'].isin(vessel_cts)]
    sc.pp.highly_variable_genes(adata_tms_vessel)
    sc.pp.pca(adata_tms_vessel, mask_var="highly_variable")
    sce.pp.harmony_integrate(adata_tms_vessel,'donor_id',max_iter_harmony = 30)
    sc.pp.neighbors(adata_tms_vessel, use_rep="X_pca_harmony")
    sc.tl.leiden(adata_tms_vessel, resolution=0.5)
    sc.tl.umap(adata_tms_vessel, min_dist=1)
    adata_tms_vessel.obs['Cell Subtype'] = pd.Categorical(adata_tms_vessel.obs['leiden'].map({'0':'Cap1','1':'Cap1','2':'Cap1','3':'Venous EC','4':'Cap1','5':'Arterial EC','6':'Cap1',}),
                                                         categories=['Arterial EC','Cap1','Venous EC'])
    adata_tms_vessel.uns['Cell Subtype_colors']= ['#4A90E2','#9B59B6','#E35D6A']

    start_and_end_states = {'root_ct': ('Cap1', 1, 'max'),
                            'terminal_cts': {'Arterial EC': (1, 'min'),
                                             'Venous EC': (0, 'min')
                                             }
                            }
    palantir_psuedotime_routine_external_datasets(adata_tms_vessel, start_and_end_states, save_prefix=figures,
                                                  adata_name=adata_name)
    find_gene_overlap_in_pseudotimes(adata_tms_vessel, save_prefix=figures, adata_name=adata_name)




'''
Goal:remove extra layers to reduce file size
conda_env:vessel_size
'''

import scanpy as sc
import scanpy.external as sce
import os 
import pandas as pd 
import numpy as np
import scipy
adata_name='venous_ec'
data = "data/single_cell_files/scanpy_files"

if __name__ == "__main__":
    adata = sc.read(f"{data}/{adata_name}_vessel_size_plot.gz.h5ad")
    ct_s_pos = adata.obs['ct_s_pos']
    del adata
    adata_all = sc.read(f'{data}/{adata_name}_celltyped_no_cc.gz.h5ad')
    adata_all.obs['Cell Subtype_w_size'] = adata_all.obs['Cell Subtype_no_cc'].astype(str).copy()
    adata_all.obs.loc[ct_s_pos.index, 'Cell Subtype_w_size'] = ct_s_pos.astype(str)
    adata_all.layers['soupx'] = adata_all.layers['raw'].copy()
    del adata_all.layers['log1p_cc_regress']
    del adata_all.layers['cp10k']
    del adata_all.uns
    del adata_all.varm
    del adata_all.obsp
    adata_all.obs = adata_all.obs[['Library', 'Treatment', 'Timepoint', 'Sort','Lineage','Cell Subtype','Cell Subtype_no_cc','Cell Subtype_w_size']]
    adata_all.var = adata_all.var[['gene_name', 'gene_id', 'gene_type', 'seqname', 'transcript_name', 'protein_id']]
    print(adata_all)
    # adata_all.write(f'{data}/vessel_size_share.gz.h5ad',compression='gzip')
    os.makedirs(os.path.join(data, "processed_txt"), exist_ok=True)
    pd.DataFrame(adata_all.var.index).to_csv(os.path.join(data, "processed_txt/genes.txt" ), sep = "\t")
    pd.DataFrame(adata_all.obs.index).to_csv(os.path.join(data, "processed_txt/barcodes.txt"), sep = "\t")
    adata_all.obs.to_csv(os.path.join(data, "processed_txt/cell_metadata.txt"), sep = "\t")
    adata_all.var.to_csv(os.path.join(data, "processed_txt/gene_metadata.txt"), sep = "\t")
    adata_all.X = adata_all.layers['raw']
    scipy.io.mmwrite(os.path.join(data, "processed_txt/matrix.mtx"), adata_all.X)





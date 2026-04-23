'''
Goal:Create figures for vessel size paper
conda_env:vessel_size
'''

import scvelo as scv
import scanpy as sc
import os
import pandas as pd 
import numpy as np
from matplotlib_venn import venn3
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy as sp
import palantir
from functions import custom_dotplot

adata_name='venous_ec'
figures = "data/figures/figures"
data = "data/single_cell_files/scanpy_files"

os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(2,2))
sc.settings.figdir = figures
sns.set_style('white', rc={
    'xtick.bottom': True,
    'ytick.left': True,
})
plt.rcParams["font.family"] = "Arial"
size=30
sc.settings.autoshow = False

if __name__ == '__main__':
    # for supplemental compare veseesl size scores by splititing treatment
    fols = {'all': 'data/figures/vessel_size_calculations',
            'Normoxia': 'data/figures/vessel_size_calculations/Normoxia',
            'Hyperoxia': 'data/figures/vessel_size_calculations/Hyperoxia',
            }
    dfs = {}
    for fol in fols:
        dfs[fol] = pd.read_csv(f'{fols[fol]}/select_obs.csv', index_col=0)

    all_df = dfs['all']
    split_df = pd.concat([dfs['Normoxia'], dfs['Hyperoxia']])

    all_df = all_df.sort_index()
    split_df = split_df.sort_index()
    delta = all_df['Vessel size score'] - split_df['Vessel size score']
    delta = pd.DataFrame(delta)
    delta['Treatment'] = ['Normoxia' if 'Nox' in x else 'Hyperoxia' for x in delta.index]

    fig, axs = plt.subplots(1, 2, figsize=(4, 2))
    axs = axs.ravel()
    plt.subplots_adjust(wspace=0.75)
    sns.regplot(x=all_df['Vessel size score'], y=split_df['Vessel size score'], scatter_kws={'s': 3, 'linewidths': 0},
                line_kws={'color': 'red'}, ci=None, ax=axs[0])
    axs[0].set_xlabel('Vessel size score\ntogether')
    axs[0].set_ylabel('Vessel size score\nsplit')
    r, p = sp.stats.pearsonr(all_df['Vessel size score'], split_df['Vessel size score'])
    axs[0].text(0.5, 1.05, 'r = {:.2f}'.format(r, p), fontsize=10,
                transform=axs[0].transAxes)
    sns.violinplot(delta, x='Treatment', y='Vessel size score', order=['Normoxia', 'Hyperoxia'],
                   palette=sns.color_palette('tab10'), ax=axs[1])
    axs[1].set_ylabel('Vessel size score delta')
    axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=45)
    axs[1].set_xlabel('')

    fig.savefig(f'{figures}/vessel_size_calculation_comparison.png', bbox_inches='tight', dpi=300)
    plt.close()
    # make venn diagram to show overlap
    fols = {'all': 'data/figures/vessel_size_calculations',
            'Normoxia': 'data/figures/vessel_size_calculations/Normoxia',
            'Hyperoxia': 'data/figures/vessel_size_calculations/Hyperoxia',
            }
    dfs = {}
    for fol in fols:
        dfs[fol] = pd.read_csv(f'{fols[fol]}/top50_shared_corr_df.csv', index_col=0)

    venn = venn3([set(dfs['all']['large_genes'].dropna().values.tolist()),
                  set(dfs['Normoxia']['large_genes'].dropna().values.tolist()),
                  set(dfs['Hyperoxia']['large_genes'].dropna().values.tolist(),
                      )],
                 set_labels=('Together', 'Normoxia', 'Hyperoxia'),
                 alpha=0.7)
    for text in venn.set_labels:
        text.set_fontsize(12)
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(12)
    plt.title("Large genes")
    plt.savefig(f'{figures}/venn_diagram_large_across_methods.png', dpi=300, bbox_inches='tight')
    plt.close()

    venn = venn3([set(dfs['all']['small_genes'].dropna().values.tolist()),
                  set(dfs['Normoxia']['small_genes'].dropna().values.tolist()),
                  set(dfs['Hyperoxia']['small_genes'].dropna().values.tolist(),
                      )],
                 set_labels=('Together', 'Normoxia', 'Hyperoxia'),
                 alpha=0.7)
    for text in venn.set_labels:
        text.set_fontsize(12)
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(12)
    plt.title("Small genes")
    plt.savefig(f'{figures}/venn_diagram_small_across_methods.png', dpi=300, bbox_inches='tight')
    plt.close()
    ## Make plots needing all cell types, mostly figure 1 and supplementals
    adata_all = sc.read(f'{data}/{adata_name}_celltyped_no_cc.gz.h5ad')
    del adata_all.layers
    ### Lineage UMAP
    adata_all.uns['Lineage_colors'] = adata_all.uns['Lineage_colors'].tolist()
    adata_all.uns['Lineage_colors'] = adata_all.uns['Lineage_colors'][::-1]
    sc.tl.umap(adata_all, min_dist=0.5)
    sc.pl.umap(adata_all,
               color='Lineage',
               alpha=0.5,
               size=5,
               frameon=False,
               title="",
               # legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=1,
               # add_outline=True,
               save='Lineage.png'
               )
    adata_all.uns['Cell Subtype_no_cc_colors'][0] = '#4A90E2' # Art
    adata_all.uns['Cell Subtype_no_cc_colors'][1] = '#9B59B6' # Cap1
    adata_all.uns['Cell Subtype_no_cc_colors'][6] = '#E35D6A' # Ven
    adata_all.uns['Cell Subtype_colors'][0] = '#4A90E2'  # Art
    adata_all.uns['Cell Subtype_colors'][1] = '#9B59B6'  # Cap1
    adata_all.uns['Cell Subtype_colors'][8] = '#E35D6A'  # Ven
    print(adata_all.uns['Cell Subtype_no_cc_colors'][:7])
    cts = ['Arterial EC','Cap1','Venous EC']
    ###Cell Subtype umap
    sc.pl.umap(adata_all,
               color='Cell Subtype_no_cc',
               groups=cts,
               na_in_legend=False,
               size=5,
               frameon=False,
               title="",
               # legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=1,
               save='celltype_highlighted.png'
               )

    adata_all.obs['Treatment'] = pd.Categorical(adata_all.obs['Treatment'],categories=['Normoxia','Hyperoxia'])
    nindex = np.random.permutation(adata_all.obs.index)
    sc.pl.umap(adata_all[nindex,:],
               color='Treatment',
               na_in_legend=False,
               size=5,
               frameon=False,
               alpha=0.7,
               title="",
               # legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=1,
               save='treatment.png'
               )
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(4,4))

    ###Cell Subtype umap
    sc.pl.umap(adata_all,
               color='Cell Subtype_no_cc',
               na_in_legend=False,
               size=5,
               frameon=False,
               title="",
               # legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=1,
               save='celltypes.png'
               )

    sc.pl.embedding(adata_all[adata_all.obs['Lineage']=='Endothelial'],
               color=['Cell Subtype'],
                    basis='X_umap_Endothelial',
               na_in_legend=False,
               size=5,
               frameon=False,
               # legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=1,
               save='endo_cts.png'
               )
    sc.pl.embedding(adata_all[adata_all.obs['Lineage']=='Endothelial'],
               color=['Cell Subtype_no_cc'],
                    basis='X_umap_Endothelial_no_cc',
               na_in_legend=False,
               size=10,
               frameon=False,
               # legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=1,
               save='endo_no_cc_cts.png'
               )
    sc.pl.embedding(adata_all[adata_all.obs['Lineage']=='Endothelial'],
               color=['Cell Subtype'],
                    basis='X_umap_Endothelial_no_cc',
                    groups = ['Proliferating Cap','Proliferating Venous EC'],
               na_in_legend=False,
               size=10,
               frameon=False,
               # legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=1,
               save='endo_no_cc_cts_pro.png'
               )
    sc.pl.DotPlot(adata_all,['Esr1','Esr2'],groupby='Cell Subtype_no_cc',).style(cmap='viridis').swap_axes().savefig(f'{figures}/dotplot_esr1_all_celltypes.png',dpi=300,bbox_inches='tight')
    ct_markers = [
    "Cdh5",
    "Gja5",
    "Kit",
    "Car4",
    "Mmrn1",
    "Col15a1",
    "Slc6a2",
    "Epcam",
    "Ager",
    "Sftpc",
    "Dnah12",
    "Scgb3a2",
    "Nrxn1",
    "Ptprc",
    "Plet1",
    "Ms4a1",
    "C1qa",
    "Mcpt8",
    "S100a4",
    "Retnlg",
    "Cd3e",
    "Itgae",
    "Mreg",
    "Col1a1",
    "Serpinf1",
    "Hhip",
    "Col13a1",
    "Dsc2",
    "Wt1",
    "Tgfbi",
    "Higd1b",
    "Nkain2",
    "Tnnt1",
    "Ctnna3"
]
    sc.pl.DotPlot(adata_all,ct_markers,standard_scale='var',groupby='Cell Subtype_no_cc',).style(cmap='viridis').swap_axes().savefig(f'{figures}/dotplot_ct_markers_all_celltypes.png',dpi=300,bbox_inches='tight')

    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
    sc.settings.figdir = figures
    sns.set_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True,
    })
    plt.rcParams["font.family"] = "Arial"
    size=30
    # all vessel size based figures
    adata = sc.read(f'{data}/{adata_name}_vessel_size_plot.gz.h5ad')
    adata.obs['Treatment'] = pd.Categorical(adata.obs['Treatment'],categories=['Normoxia','Hyperoxia'])
    sc.pl.umap(adata,
               color='Cell Subtype',
               groups=['Proliferating Cap', 'Proliferating Venous EC'],
               na_in_legend=False,
               size=size,
               frameon=False,
               title="",
               # legend_loc=None,
               legend_fontsize=10,
               legend_fontoutline=0,
               save='vascular_celltypes_pro.png'
               )
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
    ## show hyperoxia's effect on vessel size score
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
    df = sc.get.obs_df(adata, ['Treatment', 'ct_s_pos', 'Vessel size score'])
    fig, axs= plt.subplots(1, 2, figsize=(4, 2),sharey=True)
    axs=axs.ravel()
    for i,ct in enumerate(['Arterial', 'Venous']):
        if ct == 'Arterial':
            text = ['PAEC','Artery']
            df_ct = df.loc[df['ct_s_pos'].str.contains("|".join(text), case=False, na=False)]
            print(df_ct['ct_s_pos'].unique())
        else:
            text = ['PVEC', 'Vein']
            df_ct = df.loc[df['ct_s_pos'].str.contains("|".join(text), case=False, na=False)]
            print(df_ct['ct_s_pos'].unique())


        hue_order = ['Normoxia', 'Hyperoxia']
        palette = adata.uns['Treatment_colors']
        ax=axs[i]
        sns.kdeplot(
            data=df_ct,
            x="Vessel size score",
            hue='Treatment',
            hue_order=hue_order,
            palette=palette,
            common_norm=False,
            ax=ax
        )

        ax.set_title(f'{ct} branch')
        if i ==0:
            axs[i].get_legend().remove()
        else:
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(f'{figures}/histplot_vessel_size_treatment.png', dpi=300, bbox_inches='tight')
    plt.close()
    scv.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
    scv.settings.figdir = figures
    sns.set_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True,
    })
    plt.rcParams["font.family"] = "Arial"
    scv.pl.velocity_embedding_stream(adata,
                                     basis='umap',
                                     density=1,
                                     size=size,
                                     color='Cell Subtype_no_cc',
                                     legend_loc=None,
                                    title='RNA Velocity',
                                     show=False,
                                     save=f'rna_velocity.png')
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
    sc.settings.figdir = figures
    sns.set_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True,
    })
    plt.rcParams["font.family"] = "Arial"
    sc.pl.umap(adata, color='Capillary side',palette=[adata.uns['Cell Subtype_no_cc_colors'][0],adata.uns['Cell Subtype_no_cc_colors'][2]],frameon=False,alpha=0.7,size=size,na_in_legend=False,save='cap_side.png')
    fig, ax = plt.subplots(1,1,figsize=(2,2))
    palantir.plot.plot_trajectory(
        adata, # your anndata
        "Arterial EC", # the branch to plot
        cell_color="palantir_pseudotime", # the ad.obs colum to color the cells by
        n_arrows=5, # the number of arrow heads along the path
        color='#4A90E2', # the color of the path and arrow heads
        scanpy_kwargs=dict(cmap="viridis",size=size), # arguments passed to scanpy.pl.embedding
        arrowprops=dict(arrowstyle="->,head_length=.25,head_width=.25", lw=2), # appearance of the arrow heads
        lw=2, # thickness of the path
    ax=ax
        # pseudotime_interval=(0, .9), # interval of the pseudotime to cover with the path
    )
    fig.tight_layout()
    fig.savefig(f'{figures}/palantir_art_trajectory.png')
    plt.close()
    fig, ax = plt.subplots(1,1,figsize=(2,2))
    palantir.plot.plot_trajectory(
        adata, # your anndata
        "Venous EC", # the branch to plot
        cell_color="palantir_pseudotime", # the ad.obs colum to color the cells by
        n_arrows=5, # the number of arrow heads along the path
        color='#E35D6A', # the color of the path and arrow heads
        scanpy_kwargs=dict(cmap="viridis",size=size), # arguments passed to scanpy.pl.embedding
        arrowprops=dict(arrowstyle="->,head_length=.25,head_width=.25", lw=2), # appearance of the arrow heads
        lw=2, # thickness of the path
    ax=ax,
        # pseudotime_interval=(0, .9), # interval of the pseudotime to cover with the path
    )
    fig.tight_layout()
    fig.savefig(f'{figures}/palantir_ven_trajectory.png')
    plt.close()
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1, 1))
    vein_restricted = ['Slc6a2', 'Pcdh7', 'Fgl2', 'Nrg1', ]
    vein_cap = ['Nr2f2', 'Nrp2', 'Mmp16', 'Negr1']
    art_restricted = ['Gja5', 'Dkk2', 'Fmo2', 'Mannr', ]
    art_cap = ['Dll4', 'Hey1', 'Pcsk5', 'Tox2']
    sc.pl.umap(adata, color=art_restricted + vein_restricted + art_cap + vein_cap, ncols=4,
               frameon=False,
               wspace=0.2,
               hspace=0.25,
               s=size / 2.25,
               legend_loc=None,
               colorbar_loc=None,
               cmap='viridis',
               save='art_to_vein_markers_expanded.png')
    sc.set_figure_params(dpi_save=300, fontsize=8, figsize=(1,1.5))
    ## ran metascape online and got results
    pathway_szs = (('small',[2,5,6,9,14,16]),('large',[1,6,7,10,11,18]))
    for grp in pathway_szs:
        sz,ls = grp
        print(sz)
        df = pd.read_excel(f'{figures}/metascape_{sz}/metascape_result.xlsx',sheet_name='Enrichment',index_col=0,header=0)
        df = df[(df.index.str.endswith('Summary'))
           ]
        df = df.loc[df.index.to_series().str.extract(r'^(\d+)')[0].astype(int).isin(ls)] #regex to stop startswith 1 grabbing all number 10-19 etc
        df['-Log10(p-value)'] = -1 * df['LogP']
        df["Description"] = df["Description"].apply(lambda x: x[0].upper() + x[1:] if isinstance(x, str) and x else x)

        # print(df)
        fig, ax = plt.subplots(1,1,figsize=(2.5,0.75))
        ax = sns.barplot(df,x='-Log10(p-value)',
                         y='Description',
                         palette=['black']
                        )
        ax.set_ylabel('')
        sns.despine()
        fig.savefig(f'{figures}/barplot_pathway_analysis_{sz}.png',bbox_inches='tight',dpi=300)
        plt.close()
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
    sc.pl.umap(adata,color=['Vessel size score'],cmap='Oranges',frameon=False,size=size,save='vessel_size_score.png')
    sc.pl.umap(adata,color=['Vessel size category'],cmap='viridis',frameon=False,size=size,save='vessel_size_category.png')
    sc.pl.umap(adata,color=['Fbln2','Tmem100'],cmap='viridis',frameon=False,size=size,save='vessel_size_markers_validated.png')
    df = pd.read_csv('data/figures/vessel_size_calculations/top50_shared_corr_df.csv',index_col=0)
    genes = df['large_genes'].dropna().values.tolist() + df['small_genes'].dropna().values.tolist()
    fig = palantir.plot.plot_gene_trend_heatmaps(adata, genes,cmap='viridis')
    fig.tight_layout()
    fig.savefig(f'{figures}/palantir_heatmap_gene_trends.png')
    plt.close()

    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(4,4))
    sc.pl.umap(adata,color='Slc6a2',cmap='viridis',frameon=False,save='colorbar_umap.png')
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
    sns.set_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True,
    })
    plt.rcParams["font.family"] = "Arial"
    gene_ls = [
        'Ccdc85a',
        'Glp1r',
        'Kit',
             'Sox4','Nrp1','Ifitm3',
             'Ptprr','Adgrg6','Foxo1',
             'Mgp','Eln','Nr4a2',
    ]

    sc.pl.DotPlot(adata,gene_ls[::-1],groupby='ct_s_pos',categories_order=['PAEC L','PAEC M','PAEC S','Cap1-Artery', 'Cap1-Vein','PVEC S','PVEC M', 'PVEC L'],
                  standard_scale='var').style(cmap='viridis').savefig(f'{figures}/dotplot_size_markers.png',dpi=300,bbox_inches='tight')
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))

    gene_ls = ['Dkk2','Sox6',
               'Lama3','Bdkrb2',
               'Depp1','Stc1',
               'Adam23','Ntrk2',
               'Emid1','Chrm2',
               'Chrm3','Rarb',
               'Fads2b','Gria3',
               'Ptger3','Moxd1',
                 ]
    sc.pl.dotplot(adata,gene_ls,groupby='ct_s_pos',cmap='viridis',categories_order=['PAEC L','PAEC M','PAEC S','Cap1-Artery', 'Cap1-Vein','PVEC S','PVEC M', 'PVEC L'],
                  standard_scale='var',save='size_axis.png')
    sc.pl.dotplot(adata, ['Fbln2','Tmem100'], groupby='ct_s_pos', cmap='viridis',
                  categories_order=['PAEC L', 'PAEC M', 'PAEC S', 'Cap1-Artery', 'Cap1-Vein', 'PVEC S', 'PVEC M',
                                    'PVEC L'],
                  standard_scale='var', save='fbln2_tmem100.png')
    sns.set_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True,
    })

    adata.obs['Treatment_short'] = adata.obs['Treatment'].str[0]
    adata_ven = adata[adata.obs['ct_s_pos'].isin(['Cap1-Vein','PVEC L', 'PVEC M', 'PVEC S'])].copy()
    df = sc.get.obs_df(adata_ven,['proliferation_score','Vessel size category', 'Treatment'])
    df.rename(columns={'proliferation_score':'Proliferation score'},inplace=True)
    fig, axs = plt.subplots(1, 2,figsize=(3,2),sharey=True)
    axs = axs.ravel()
    hue_order = ['capillary','small', 'medium', 'large']
    palette = adata.uns['Vessel size category_colors']

    for i, treat in enumerate(['Normoxia', 'Hyperoxia']):
        ax = sns.kdeplot(
            data=df[df['Treatment'] == treat],
            x="Proliferation score",
            hue='Vessel size category',
            hue_order=hue_order,
            palette=palette,
            fill=False,
            common_norm=False,
            ax=axs[i]
        )

        ax.set_title(treat)
        ax.set_ylabel('Density')
        ax.set_xlabel('')
        ax.set_xticks([0, 0.75])
        ax.set_xticklabels(['0', '0.75'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if i == 0:
            handles = [
                mlines.Line2D([], [], color=palette[j], label=label, linewidth=2)
                for j, label in reversed(list(enumerate(hue_order)))
            ]
            ax.legend(handles=handles[::-1], frameon=False, fontsize="8", title='', loc="upper right")
        else:
            ax.get_legend().remove()

    fig.supxlabel('Proliferation score', y=0.15, x=0.52)
    fig.tight_layout()
    fig.savefig(f'{figures}/histplot_venous_ec_treat_proliferation_score.png', dpi=300, bbox_inches='tight')
    plt.close()
    adata_ven.obs['Treatment'] = adata_ven.obs['Treatment'].str[0]
    adata_art = adata[adata.obs['ct_s_pos'].isin(['Cap1-Artery','PAEC L', 'PAEC M', 'PAEC S'])].copy()
    adata_art.obs['Treatment'] = adata_art.obs['Treatment'].str[0]

    sns.set_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True,
    })
    custom_dotplot(adata_ven,
                   ['Mki67', 'Top2a', 'Birc5', 'Hmgb2', 'Cenpf'],
                   scale_by_gene=True,
                   x_obs='Treatment_short',
                   y_obs='Vessel size category',
                   x_order=['N','H'],
                   cmap='viridis',
                   dot_max_size=100,
                   figsize=(4,2),
                   save=f'{figures}/dotplot_proliferation_vein.png')

    for gene_ls in [
    ["Thsd7a", "Sox7", "Cd44", "Igf1r"],
    ["Acvr1b", "Inhbb", "Smurf1", "Ccm2"],
    ["Serpine2", "Pros1", "Sdc4", "Bax", "Mdm2"],
]:
        custom_dotplot(
            adata[adata.obs['ct_s_pos'].isin(['Cap1-Artery', 'Cap1-Vein'])].copy(),
            gene_ls,
            scale_by_gene=True,
            x_obs='Treatment_short',
            y_obs='ct_s_pos',
            y_order=['Cap1-Artery', 'Cap1-Vein'],
            x_order=['N', 'H'],
            cmap='viridis',
            dot_max_size=100,
            figsize=(4.5, 1.5),
            save=f'{figures}/dotplot_custom_hyperoxia_genes_caps_{gene_ls[0]}.png')

    sc.pl.DotPlot(adata,['Mki67', 'Top2a', 'Birc5', 'Hmgb2', 'Cenpf'],groupby='ct_s_pos',cmap='viridis',categories_order=['PAEC L','PAEC M','PAEC S','Cap1-Artery', 'Cap1-Vein','PVEC S','PVEC M', 'PVEC L'],
                  standard_scale='var').savefig(f'{figures}/dotplot_size_axis_proliferation_all.png',bbox_inches='tight',dpi=300)
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1,1))
    sc.settings.figdir = figures
    sns.set_style('white', rc={
        'xtick.bottom': True,
        'ytick.left': True,
    })
    size=15
    plt.rcParams["font.family"] = "Arial"
    sc.pl.umap(adata,color=['Cxcl12','Cxcr4','Ackr3'],colorbar_loc=None,wspace=0.2,hspace=0.25,ncols=2,s=size/4,frameon=False,cmap='viridis',save='cxcl12_signaling.png')
    sc.pl.umap(adata,color=['Cxcl12','Cxcr4','Ackr3'],wspace=0.2,hspace=0.25,ncols=2,s=size/2.25,frameon=False,cmap='viridis',save='cxcl12_signaling_cb.png')
    sc.pl.DotPlot(adata[adata.obs['ct_s_pos'].isin(['PAEC L','PAEC M','PAEC S','Cap1-Artery',])],
                  ['Cxcl12','Cxcr4'],
                  standard_scale='var',
                  groupby='Vessel size category').style(cmap='viridis').savefig(f'{figures}/dotplot_art_cxcl12_signal.png',dpi=300)
    sc.pl.DotPlot(adata[adata.obs['ct_s_pos'].isin(['PVEC L','PVEC M','PVEC S','Cap1-Vein',])]
                  ,['Ackr3'],
                  standard_scale='var',groupby='Vessel size category').style(cmap='viridis').savefig(f'{figures}/dotplot_ven_cxcl12_signal.png',dpi=300)
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
    sc.pl.umap(adata,color='Esr2',frameon=False,s=size,cmap='viridis',save='_Esr2.png')
    sc.set_figure_params(dpi_save=300, fontsize=10, figsize=(1.5,1.5))
    sc.pl.DotPlot(adata,
                  ['Esr1'],
                  groupby = 'ct_s_pos',
                  categories_order=['PAEC L','PAEC M','PAEC S','Cap1-Artery', 'Cap1-Vein','PVEC S','PVEC M', 'PVEC L']
                 ).style(cmap='viridis').savefig(f'{figures}/dotplot_ec_esr1.png',dpi=300)

    # CCC plots
    df_ccc=pd.read_csv(f'{figures}/liana_vessel_size_ccc.csv')

    artery_ligands = ['Dkk2', 'Efna5', 'Lama3', 'Jag2']
    artery_side = ['PAEC L', 'PAEC M', 'PAEC S', 'Cap1-Artery']
    vein_ligands = ['Col18a1','Vcan','Jam2','Dll1']
    vein_side = ['PVEC L','PVEC M','PVEC S', 'Cap1-Vein',]
    caps = [ 'Cap1-Artery', 'Cap1-Vein']
    cap_ligands = ['Adam23']
    cap_partners = ['Cap1-Artery','Cap1-Vein','Cap2',
               'Pericyte']
    mural = ['Vascular smooth muscle', 'Pericyte']
    for x in ['PAEC','PVEC','Cap']:
        if x == 'PAEC':
            side = artery_side
            ligs = artery_ligands
        elif x == 'PVEC':
            side = vein_side
            ligs = vein_ligands
        elif x == 'Cap':
            side = caps
            ligs = cap_ligands
        sc.pl.DotPlot(adata[adata.obs['ct_s_pos'].isin(side)],
                      ligs,
                      groupby='ct_s_pos',
                      categories_order=side,
                      standard_scale='var'
                      ).style(cmap='viridis',largest_dot=150,dot_max=1).savefig(f'{figures}/dotplot_{x}_ligands_by_size.png', dpi=300)
        for lig in ligs:
            complexes = df_ccc.loc[(df_ccc.source.isin(side)) &
                                       (df_ccc.target.isin(mural)) &
                                       (df_ccc.ligand_complex.isin([lig]))
                                       ]['receptor_complex'].unique().tolist()
            receptor_genes = []
            for receptor in complexes:
                if "_" in receptor:
                    receptor_genes.extend(receptor.split("_"))
                else:
                    receptor_genes.append(receptor)
            receptor_genes = sorted(set(receptor_genes))

            adata_all.obs['Cell Subtype_w_size'] = adata_all.obs['Cell Subtype_no_cc'].astype(str).copy()
            adata_all.obs.loc[adata.obs.index, 'Cell Subtype_w_size'] = adata.obs['ct_s_pos'].astype(str)
            if x in ['PAEC','PVEC']:
                sc.pl.DotPlot(adata_all[adata_all.obs['Cell Subtype_no_cc'].isin(mural)],
                              receptor_genes,
                              groupby='Cell Subtype_no_cc',
                              title=f'Receptors for {lig} in mural',
                              ).style(cmap='viridis',largest_dot=150,dot_max=1).savefig(f'{figures}/dotplot_{x}_{lig}_mural_receptors.png', dpi=300)
            elif x in ['Cap']:
                sc.pl.DotPlot(adata_all[adata_all.obs['Cell Subtype_w_size'].isin(cap_partners)],
                              receptor_genes+['Itgav'],
                              groupby='Cell Subtype_w_size',
                              title=f'Receptors for {lig} in alveoli vasculature',
                              ).style(cmap='viridis', largest_dot=150, dot_max=1).savefig(
                    f'{figures}/dotplot_{x}_{lig}_alv_vasc_receptors.png', dpi=300)
        sc.pl.DotPlot(adata_all[adata_all.obs['Cell Subtype_no_cc'].isin(mural)],
                      ['Epha1','Epha3','Epha4'],
                      groupby='Cell Subtype_no_cc',
                      title=f'Receptors for Efna5 in mural',
                      ).style(cmap='viridis', largest_dot=150, dot_max=1).savefig(
            f'{figures}/dotplot_Efna5_mural_receptors_custom.png', dpi=300)
        sc.pl.DotPlot(adata_all[adata_all.obs['Cell Subtype_no_cc'].isin(mural)],
                      ['Itga3','Itga6','Itgb1','Sdc2'],
                      groupby='Cell Subtype_no_cc',
                      title=f'Receptors for Lama3 in mural',
                      ).style(cmap='viridis', largest_dot=150, dot_max=1).savefig(
            f'{figures}/dotplot_Lama3_mural_receptors_custom.png', dpi=300)
        sc.pl.DotPlot(adata_all[adata_all.obs['Cell Subtype_no_cc'].isin(mural)],
                      ['Itga3','Itga4','Itgav','Itgb1'],
                      groupby='Cell Subtype_no_cc',
                      title=f'Receptors for  Col18a1/Jam2 in mural',
                      ).style(cmap='viridis', largest_dot=150, dot_max=1).savefig(
            f'{figures}/dotplot_Col18a1_Jam2_mural_receptors_custom.png', dpi=300)


        def annotate(data, x, y, ax, x_loc, y_loc, **kws):
            r, p = sp.stats.pearsonr(np.log10(data[x]), data[y])
            # r*=r
            ax.text(x_loc, y_loc, 'r = {:.2f}\np = {:.2g}'.format(r, p), fontsize=10,
                    transform=ax.transAxes)


        sns.set_style('white', rc={
            'xtick.bottom': True,
            'ytick.left': True,
        })
        plt.rcParams["font.family"] = "Arial"
        ## EdU images
        data = '/media/carsten/hdd/documents/manuscripts/vessel_size/images/TIF proliferation/intensity_output/Image.csv'
        df = pd.read_csv(data, header=0)
        df = df.loc[df['Intensity_TotalIntensity_DAPI'] > 100]  # Takes out misclipped vessels
        df['Replicate'] = df['FileName_composite'].str[0]
        df['Treatment'] = df['FileName_composite'].str[1]
        df['Treatment'] = df['Treatment'].replace(' ', 'N')
        df['Mouse'] = df['Treatment'] + df['Replicate']
        df['EdU/DAPI (AU)'] = df['Intensity_MeanIntensity_EdU'] / df['Intensity_MeanIntensity_DAPI']
        df['Slc6a2/DAPI (AU)'] = df['Intensity_MeanIntensity_SLC6A2'] / df['Intensity_MeanIntensity_DAPI']
        df['EdU/SLC6A2 (AU)'] = df['Intensity_MeanIntensity_EdU'] / df['Intensity_MeanIntensity_SLC6A2']
        df['% EdU+ PVEC'] = df['Count_EdU_Slc6A2'] / df['Count_nuc_slc6a2'] * 100
        df['Diameter (µm)'] = df[['Height_DAPI', 'Width_DAPI']].T.max() * .137
        df = df.dropna(subset=['Diameter (µm)', '% EdU+ PVEC', 'EdU/DAPI (AU)'])
        df['diameter_bin'] = pd.cut(df['Diameter (µm)'], bins=[0, 50, 1000], labels=['≤50', '>50'])
        df2 = df.groupby(['Mouse', 'diameter_bin'])[['Count_EdU_Slc6A2', 'Count_nuc_slc6a2']].sum()
        df2['Treatment'] = df2.index.get_level_values(0).str[0]
        df2['% EdU+ PVEC'] = df2['Count_EdU_Slc6A2'] / df2['Count_nuc_slc6a2'] * 100
        df2.rename_axis(index=['Mouse', 'Diameter (µm)'], inplace=True)
        df2.reset_index(inplace=True)
        fig, axs = plt.subplots(1, 2, figsize=(3.5, 1.5))
        plt.subplots_adjust(wspace=0.35)
        axs = axs.ravel()
        sns.regplot(df, x=df['Diameter (µm)'], y=df['% EdU+ PVEC'], scatter_kws={'s': 10, 'linewidths': 0},
                    line_kws={'color': 'red'}, logx=True, ci=None, ax=axs[0])
        for i in [0]:
            axs[i].set_xscale('log')
            axs[i].set_xticklabels(['', '', 100])
            axs[i].set_ylim([0, axs[i].get_ylim()[1]])

        annotate(df, 'Diameter (µm)', '% EdU+ PVEC', axs[0], 0.5, 1.05)
        sns.barplot(df2, x='Diameter (µm)', y='% EdU+ PVEC', hue='Treatment', ax=axs[1], hue_order=['N', 'H'], alpha=0.3,
                    ci=None)
        sns.stripplot(df2, x='Diameter (µm)', y='% EdU+ PVEC', hue='Treatment', size=5, dodge=True, ax=axs[1],
                      hue_order=['N', 'H'], edgecolor='black', linewidth=0.5)
        for ax in axs:
            ax.spines[['right', 'top']].set_visible(False)

        handles, labels = axs[1].get_legend_handles_labels()
        unique_labels = list(dict.fromkeys(labels))
        unique_handles = [handles[labels.index(lbl) + 2] for lbl in unique_labels]
        axs[1].legend(unique_handles, unique_labels, title='Treatment', loc='upper left', bbox_to_anchor=(1, 1))
        axs[1].legend_.set_frame_on(False)

        i = 0
        for size in df2['Diameter (µm)'].unique():
            df_test = df2.loc[df2['Diameter (µm)'] == size]
            norm = df_test[df_test['Treatment'] == 'N']['% EdU+ PVEC']
            hyper = df_test[df_test['Treatment'] == 'H']['% EdU+ PVEC']

            t_statistic, p_value = sp.stats.ttest_ind(norm, hyper, nan_policy='omit', equal_var=True)
            print(f'In {size}')
            print(f"Independent t-test:")
            print(f"T-statistic: {t_statistic}")
            print(f"P-value: {p_value}")
            y_coord = df2['% EdU+ PVEC'].max() + 0.25  # Adjust this value as needed
            axs[1].text(x=i, y=y_coord, s=f'p\n{p_value:.3f}', ha='center', va='bottom', fontsize=10, color='black')
            i += 1
        axs[1].set_ylabel('')
        fig.savefig(f'{figures}/scatterplots_vein_EdU.png', dpi=300, bbox_inches='tight')
        fig.show()

        ## Artery images
        data = '/media/carsten/hdd/documents/manuscripts/vessel_size/images/TIF artery/intensity_output/Image.csv'
        df = pd.read_csv(data, header=0)
        df = df.loc[df['Intensity_TotalIntensity_DAPI'] > 100]  # Takes out misclipped vessels
        df['Mouse'] = df['FileName_composite'].str[0]
        df['Dkk2/DAPI (AU)'] = df['Intensity_MeanIntensity_large_marker'] / df['Intensity_MeanIntensity_DAPI']
        df['Gja5/DAPI (AU)'] = df['Intensity_MeanIntensity_vessel_marker'] / df['Intensity_MeanIntensity_DAPI']
        df['Diameter (µm)'] = df[['Height_vessel_marker', 'Width_vessel_marker']].T.max() * .137
        fig, axs = plt.subplots(1, 2, figsize=(3.5, 1.5))
        plt.subplots_adjust(wspace=1)
        axs = axs.ravel()
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Dkk2/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[0])
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Gja5/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[1])
        for i in [0, 1,]:
            axs[i].set_xscale('log')
            axs[i].set_xticklabels(['', '', 100])
            axs[i].set_ylim([0, axs[i].get_ylim()[1]])

        annotate(df, 'Diameter (µm)', 'Dkk2/DAPI (AU)', axs[0], 0.5, 1.05)
        annotate(df, 'Diameter (µm)', 'Gja5/DAPI (AU)', axs[1], 0.5, 1.05)
        for ax in axs:
            ax.spines[['right', 'top']].set_visible(False)
        fig.savefig(f'{figures}/scatterplot_artery.png', dpi=300, bbox_inches='tight')
        fig.show()

        ## vein images
        data = '/media/carsten/hdd/documents/manuscripts/vessel_size/images/TIF vein/intensity_output/Image.csv'
        df = pd.read_csv(data, header=0)
        df = df.loc[df['Intensity_TotalIntensity_DAPI'] > 100]  # Takes out misclipped vessels
        df['Mouse'] = df['FileName_composite'].str[0]
        df['Moxd1/DAPI (AU)'] = df['Intensity_MeanIntensity_large_marker'] / df['Intensity_MeanIntensity_DAPI']
        df['Slc6a2/DAPI (AU)'] = df['Intensity_MeanIntensity_vessel_marker'] / df['Intensity_MeanIntensity_DAPI']
        df['Diameter (µm)'] = df[['Height_vessel_marker', 'Width_vessel_marker']].T.max() * .137

        fig, axs = plt.subplots(1, 2, figsize=(3.5, 1.5))
        plt.subplots_adjust(wspace=1)
        axs = axs.ravel()
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Moxd1/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[0])
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Slc6a2/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[1])

        for i in [0, 1,]:
            axs[i].set_xscale('log')
            axs[i].set_xticklabels(['', '', 100])
            axs[i].set_ylim([0, axs[i].get_ylim()[1]])

        annotate(df, 'Diameter (µm)', 'Moxd1/DAPI (AU)', axs[0], 0.5, 1.05)
        annotate(df, 'Diameter (µm)', 'Slc6a2/DAPI (AU)', axs[1], 0.5, 1.05)
        for ax in axs:
            ax.spines[['right', 'top']].set_visible(False)
        fig.savefig(f'{figures}/scatterplot_vein.png', dpi=300, bbox_inches='tight')
        fig.show()

        # all vessels
        data = '/media/carsten/hdd/documents/manuscripts/vessel_size/images/tmem_fbln2_cellprofiler/cellprofiler_output_spreadsheets_1/Image.csv'
        df = pd.read_csv(data, header=0)
        df = df.loc[df['Intensity_TotalIntensity_DAPI'] > 100]  # Takes out misclipped vessels
        df['Mouse'] = df['FileName_DAPI'].str[0]
        df['Fbln2/DAPI (AU)'] = df['Intensity_MeanIntensity_Fbln2'] / df['Intensity_MeanIntensity_DAPI']
        df['Tmem100/DAPI (AU)'] = df['Intensity_MeanIntensity_Tmem100'] / df['Intensity_MeanIntensity_DAPI']
        df['Tmem100/DAPI (AU)'] = df['Intensity_MeanIntensity_Tmem100'] / df['Intensity_MeanIntensity_DAPI']
        df['Pecam1/DAPI (AU)'] = df['Intensity_MeanIntensity_Cd31'] / df['Intensity_MeanIntensity_DAPI']
        df['Fbln2/Tmem100 (AU)'] = np.log10(df['Intensity_MeanIntensity_Fbln2'] / df['Intensity_MeanIntensity_Tmem100'])
        df['Diameter (µm)'] = df[['Height_Cd31', 'Width_Cd31']].T.max() * .137
        df = df.loc[df['Mouse'] != '4']  # mouse 4 has low quality in Fbln2 channel
        blacklist = ['3NP3_Cd31FITCFbln2Cy3Tmem100Cy5_20X_-13Apo_EDFMIP-30-Image Export-30_c2_vessel_all_3.tiff',
                     '3NP3_Cd31FITCFbln2Cy3Tmem100Cy5_20X_-13Apo_EDFMIP-30-Image Export-01_c2_vessel_all_4.tiff'
                     ]  # miscropped images

        df = df.loc[~df['FileName_Cd31'].isin(blacklist)]
        fig, axs = plt.subplots(1, 3, figsize=(5.25, 1.5))
        plt.subplots_adjust(wspace=1)
        axs = axs.ravel()
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Tmem100/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[0])
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Fbln2/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[1])
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Pecam1/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[2])

        for i in [0, 1,2]:
            axs[i].set_xscale('log')
            axs[i].set_xticklabels(['', '', 100])
            axs[i].set_ylim([0, axs[i].get_ylim()[1]])

        annotate(df, 'Diameter (µm)', 'Tmem100/DAPI (AU)', axs[0], 0.5, 1.05)
        annotate(df, 'Diameter (µm)', 'Fbln2/DAPI (AU)', axs[1], 0.5, 1.05)
        annotate(df, 'Diameter (µm)', 'Pecam1/DAPI (AU)', axs[2], 0.5, 1.05)
        for ax in axs:
            ax.spines[['right', 'top']].set_visible(False)
        fig.savefig(f'{figures}/scatterplots_fbln2_tmem100_diameter.png', dpi=300, bbox_inches='tight')
        fig.show()

        ## Adgrg6 images
        data = '/media/carsten/hdd/documents/manuscripts/vessel_size/images/Adgrg6 TIF/intensity_output/Image.csv'
        df = pd.read_csv(data, header=0)
        df = df.loc[df['Intensity_TotalIntensity_DAPI'] > 100]  # Takes out misclipped vessels
        df['Mouse'] = df['FileName_composite'].str[0]
        df['Adgrg6/DAPI (AU)'] = df['Intensity_MeanIntensity_medium_marker'] / df['Intensity_MeanIntensity_DAPI']
        df['Pecam1/DAPI (AU)'] = df['Intensity_MeanIntensity_vessel_marker'] / df['Intensity_MeanIntensity_DAPI']
        df['Adgrg6/Pecam1 (AU)'] = df['Intensity_MeanIntensity_medium_marker'] / df['Intensity_MeanIntensity_vessel_marker']

        df['Diameter (µm)'] = df[['Height_vessel_marker', 'Width_vessel_marker']].T.max() * .137
        fig, axs = plt.subplots(1, 2, figsize=(3.5, 1.5))
        plt.subplots_adjust(wspace=1)
        axs = axs.ravel()
        sns.scatterplot(df, x=df['Diameter (µm)'], y=df['Adgrg6/DAPI (AU)'],s=10, linewidths = 0,  ax=axs[0])
        sns.regplot(df, x=df['Diameter (µm)'], y=df['Pecam1/DAPI (AU)'], logx=True, ci=None,
                    scatter_kws={'s': 10, 'linewidths': 0}, line_kws={'color': 'red'}, ax=axs[1])
        for i in [0, 1,]:
            axs[i].set_xscale('log')
            axs[i].set_xticklabels(['', '', 100])
            axs[i].set_ylim([0, axs[i].get_ylim()[1]])

        annotate(df, 'Diameter (µm)', 'Pecam1/DAPI (AU)', axs[1], 0.5, 1.05)

        for ax in axs:
            ax.spines[['right', 'top']].set_visible(False)
        fig.savefig(f'{figures}/scatterplot_adgrg6_images.png', dpi=300, bbox_inches='tight')
        fig.show()

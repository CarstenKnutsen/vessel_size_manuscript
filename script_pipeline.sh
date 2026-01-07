#!/bin/bash
echo 'starting'
source /home/carsten/anaconda3/etc/profile.d/conda.sh
conda activate soupxR
Rscript soupX.R
conda activate venous_ec
python initial_analysis.py
python cell_typing_by_lineage.py
python cell_typing_by_lineage_no_cc.py
conda activate vessel_size
python create_count_table_velocyto.py
python vessel_size_trajectory.py
python vessel_size_calculations.py
conda activate liana
python vessel_size_ccc.py
conda activate vessel_size
python vessel_size_figures.py
python export_clean_h5ad.py
python hurskainen_vessel_size.py
python tabula_muris_senis_vessel_size.py
python bhattacharya_vessel_size.py
python lungmap_vessel_size.py
python tabula_sapiens_vessel_size.py
python integrate_all_vessel_size.py




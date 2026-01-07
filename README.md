# Code used in Sveiven and Knutsen et al. on single cell RNA sequencing of endothelial zonation along the pulmonary vascular tree

Alignment was done on Stanford's HPC Sherlock
- alignment.sbatch
  - Align libraries using cellranger

Conda environment files can be found in conda_envs
- functions.py
  - contains functions used throughout repo
  
All subsquent analysis follows in the order of script_pipeline.sh

- soupX.R
  - Uses SoupX to clean up ambient RNA from alignment done by cellranger
- initial_analysis.py
  - Reading in of counts, QC steps and categorizing lineages
- cell_typing_by_lineage.py
  - Reclustering and embedding lineages to assign cell types
- cell_typing_by_lineage_no_cc.py
  - Repeat of above with cell cycle genes regressed
- create_count_table_velocyto.py
  - Create count tables for all cells from velocyto alignment for RNA velocity
- vessel_size_trajectory.py
  - Run RNA velocity and pseudotime across vascular cell types
- vessel_size_calculations.py
  - Calculate estimated vessel size score
- vessel_size_ccc.py
  - Use liana to infer cell-cell communication using vessel sized cell types
- vessel_size_figures.py
  - Create figures used in Sveiven and Knutsen et al.
- export_clean_h5ad.py
  - Export files for sharing
- External datasets
  - Files below are to run each individual dataset through vessel size scoring
    - hurskainen_vessel_size.py
    - tabula_muris_senis_vessel_size.py
    - bhattacharya_vessel_size.py
    - lungmap_vessel_size.py
    - tabula_sapiens_vessel_size.py
- integrate_all_vessel_size.py
  - Integrate all datasets vessel size scoring, create figures for Sveiven and Knutsen et al.





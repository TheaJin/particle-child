Description:
   This script processes particle deposition data from simulation experiments
   on children (virtual and constricted conditions) and clustered data. It:
     • Reads and processes virtual child data (with various constriction levels)
      and cluster data.
     • Performs paired comparisons (e.g., paired t-test of intra-thoracic deposition
#       between VC (1.0) and Cluster_12).
#     • Generates a series of ggplot2 visualizations for different deposition measures:
#         - Intra-thoracic (TDF)
#         - Bronchial (BDE)
#         - Alveolar (ADF)
#     • Merges adult simulation data (from TDF, BDF, ADF sheets) with virtual child data
#       for multi-measure comparisons.
#     • Runs pairwise comparisons using emmeans and computes compact letter displays (CLD).
#     • Optionally, reads external "lobe summary" files, summarizes by condition and lobe,
#       and plots deposition by lobe.
#
# Data Requirements:
#   - Excel files containing virtual child deposition data (e.g., "particle_deposition_results_virtual_old.xlsx",
#     "particle_deposition_results_constrict_90.xlsx", "particle_deposition_results_constrict_85.xlsx",
#     "particle_deposition_results_constrict_80.xlsx").
#   - An Excel file for cluster deposition data (e.g., "particle_deposition_results_cluster_12.xlsx").
#   - An adult simulation Excel file ("result_new.xlsx") with sheets "TDF_Lung", "BDF", "ADF".
#   - A demographic dataset ("demo") with subject IDs, gender, etc.

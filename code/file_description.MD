**calculate_tpr_predict_genes.py** <- This is the script run after 'functionalML.R' script to get final predicted files. It takes in the 10 input 
files 'pred_scores*.csv'. And gives output file 'predicted_genes.csv'.

**clustering_functions.R** <- This is the script run to get the cluster of functions. It takes in the 2 input 
files 'tf_mat_tpr_clustering.csv' and 'tf_cor_tpr_clustering.csv'. And gives output file 'cluster_of_functions.pdf'.

**functionalML.R** <- The script is run to get prediction files for 295 gene-sets. It takes in 3 types of input files,
10 files of '/inst/extdata/peakall*.rds', remaining 'pathway-genes.gmt.csv' and  'unionpeak2.0'.
It gives 6 output files: 'f-1-295.txt', 'TPR-1-295.txt', 'FDR-1-295.txt', 'pred-1-295.txt', 'tfs-1-295.txt' and 'cf-1-295.txt'.

**lncrna_crispr_validation.R** <- This script reproduces results of the paper. It takes in 2 input files:
'crispr_lncrna_cell_cycle_gm128.csv', 'crispr_lncrna_cell_cycle_k562.csv'. 
It gives 2 output files: 'compare_cluster_random_lncrna_gm12878.pdf' and 'compare_cluster_random_lncrna_k562.pdf'.

**package_test_crispr.R** <- This cript is to test the R package 'GFPredict'. 
It takes in 6 input files:
'cell_cycle_screen.csv', 'crispr_lncrna_cell_cycle_gm128.csv', 'crispr_lncrna_cell_cycle_k562.csv', 'peptide_accu_screen.csv',
'resis_bacteria_screen.csv', 'viability_screen.csv'.
It gives out 6 output files:
'compare_model_null_cell_cycle.pdf', 'compare_model_null_cell_cycle_lncrna_gm12878', 'compare_model_null_cell_cycle_lncrna_k562.csv',
'compare_model_null_peptide_accu.pdf', 'compare_model_null_resis_bacteria.pdf', 'compare_model_null_viability.pdf'.

**plot_abtract_mining_results.R** <- This is run after 'pubmed_abtract_mining.py'. This will output barplot to show abtract mining results.
It takes in 2 input files: 'term_freq.csv', 'null_freq.csv'. It gives out 1 output file, 'box_plot_cooccurrence.pdf'.


**pubmed_abtract_mining.py** <- This is run to get the abstract mining of predicted gene terms and function terms co-occrrence count.
Takes in 3 files as input 'pathway-genes.gmt.csv', 'unionpeak2.0', especially 'predicted_genes.csv' file.
It gives out 4 files: 'term_freq.csv', 'pmid_mat.csv', 'null_freq.csv', 'null_pmid.csv'. 

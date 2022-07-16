Files are orderd with respect to their usage in the different scripts.

**pred_scores*.csv (10 files)** <- Input of calculate_tpr_predict_genes, contains probablity scores for all the genes (rows). 
**predicted_genes.csv** <- Output of calculate_tpr_predict_genes, contains contains predicted genes (columns) for each functions(rows). 

---


**tf_mat_tpr_clustering.csv** <-  This is the input for 'clustering_functions.R', 
contains top predictors of random forest model for each funciton.

**tf_cor_tpr_clustering.csv** <- This is the input for 'clustering_functions.R', 
contains correlation scores for top predictors of random forest model for each funciton.

**cluster_of_functions.pdf** <- This is output file of 'clustering_functions.R'. It is the cluster plot.

---

**pathway-genes.gmt.csv** <- Gene-set file, input of 'functionalML.R' and 'pubmed_abtract_mining.py'
**unionpeak2.0** <- Curated list of promoters. Input of 'functionalML.R' and 'pubmed_abtract_mining.py'

---

**f-1-295.txt** <- Output file of functionalML.R, it is the correlation value between result and ground truth used in the model taining for each function (rows).
**TPR-1-295.txt** <- Output file of functionalML.R. True positive rate for each 292 functions (columns).
**FDR-1-295.txt** <- Output file of functionalML.R. False positive rate for each 292 functions (columns).
**pred-1-295.txt** <- Output file of functionalML.R and input file for 'calculate_tpr_predict_genes.py'. Probability scores for each gene (rows), for 295 functions (columns).
**tfs-1-295.txt** <- Output file of functionalML.R. Top predictors for each function (columns).
**cf-1-295.txt** <- Output file of functionalML.R. Model performance metric for each function (rows). 

---

**crispr_lncrna_cell_cycle_gm128.csv** <- Input file for 'lncrna_crispr_validation.R' and 'package_test_crispr.R'. It is CRISPR screen for lncRNA cell cycle.
**crispr_lncrna_cell_cycle_k562.csv** <- Input file for 'lncrna_crispr_validation.R' and 'package_test_crispr.R'. It is CRISPR screen for lncRNA cell cycle.

---

**cell_cycle_screen.csv** <- Input file for 'package_test_crispr.R', it is CRISPR screen for cell cycle.
**peptide_accu_screen.csv** <- Input file for 'package_test_crispr.R', it is CRISPR screen for peptide accumulation.
**resis_bacteria_screen.csv** <- Input file for 'package_test_crispr.R', it is CRISPR screen for resistance to bacteria.
**viability_screen.csv** <- Input file for 'package_test_crispr.R', it is CRISPR screen for viability.

---

**term_freq.csv** <- This is input file for plot_abtract_mining_results.R. It is count of co-occurrence of predicted genes and corresponding function term for 295 functions. Out of which 49 functions had co-occrrence.
**null_freq.csv** <- This is input file for plot_abtract_mining_results.R. It is count of co-occurrence of predicted genes and corresponding function term for 49 functions. Out of which 7 functions had co-occrrence.

---

**predicted_genes.csv** <- Input file for 'pubmed_abtract_mining.py' and output file of 'calculate_tpr_predict_genes.py'. It contins predicted genes for each function (rows).

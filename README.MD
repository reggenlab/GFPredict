## Introduction[![](./extras/images/pin.svg)](#introduction)

**GFPredict** is a R package. It predicts coding and non-coding genes related to (user input) gene-sets. GFPredict compares the similarity between the epigenomic and transcription factor binding patterns at the promoter sites of the input genes and other genes, to make predictions. [Click here](https://www.biorxiv.org/content/10.1101/2022.04.17.488570v1) to learn more about GFPredict.

---

## Table of contents[![](./extras/images/pin.svg)](#table-of-contents)
1. [GFPredict Package installation](#installation)

2. [GFPredict usage examples](#examples)

3. [Generate results used in the analyses of the paper](#results)
   - [Results of random forest algorithm for 5917 ontological gene-sets](#ontology_results)
   - [PUBMED-based abstract mining](#pubmed_abtract_mining)
   - [Clustering of functions](#function_cluster)
   - [Validation of predicted lncRNAs in clusters using CRISPR screens](#cluster_crispr_screen)
   - [GFPredict testing using different CRISPR screens](#different_crispr_screen)


---

## 1. GFPredict Package installation[![](./extras/images/pin.svg)]<a name="installation"></a>

### Download zip file and install the package:

**Step 1.** Download zip file.

```shell
wget https://github.com/reggenlab/GFPredict/archive/refs/heads/main.zip

```
**Step 2.** Install the package in R session.

```r
devtools::install_local('main.zip')

```

### Or use following installation through R commands.


``` r
options(timeout=9999999)
devtools::install_github("reggenlab/GFPredict")
```

In case of difficulty in installation of libraries,"randomForest", "xgboost" or "glmnet", run the following commands:
``` r
packageurl_xgboost <- "http://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_0.90.0.2.tar.gz"
install.packages(packageurl_xgboost, repos=NULL, type="source") #for xgboost

install_version("glmnet", version = "3.0-2") #https://cran.r-project.org/src/contrib/Archive/glmnet/ 

packageurl_rf <- "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz"
install.packages(packageurl_rf, repos=NULL, type="source") 

```

---

## 2. GFPredict usage examples [![](./extras/images/pin.svg)]<a name="examples"></a>

After the installation:

``` r
library(GFPredict)

cell_cycle_related_genes <- c("RPL23A", "IPO4", "TARS1", "COG4", "TEX10", "TRMT112", "TXNL4A", "CLP1", 
"HSPA5", "ANAPC4", "RNF168", "ATP5F1D", "RUVBL2", "NMT1", "PNKP", "SF3B3", "FDPS", "FARSB", "HARS1",
"RPL8", "CCT3", "AQR", "MCL1", "CENPM") 

result <- predict_related_genes(cell_cycle_related_genes) #Default parameter are selected

#Fine tune the model by changing the number of bootstraps of the negative set
result <- predict_related_genes(genes = cell_cycle_related_genes, n_bootstrap = 5)

#Use different ML algorithm if the model performance is not good after using the default random forest algoritm 
result <- predict_related_genes(genes = cell_cycle_related_genes, ml_model = 'svm', n_bootstrap = 10, feature_type = 'transcription_factors')

result[[1]] # contains list of predicted genes 

result[[2]] # contains table of model performance metrics 

result[[3]] # contains table of top predictors (if ml_model = random.forest)

```

---

## 3. Generate results used in the analyses of the paper [![](./extras/images/pin.svg)] <a name="results"></a>

Run the following scrips to reproduce the results showed in the [paper](https://www.biorxiv.org/content/10.1101/2022.04.17.488570v1)

### Results of random forest algorithm for 295 ontological gene-sets <a name="ontology_results"></a>

**Step 1.**

Download the following files and set their location path inside ./code/functionML.R script:

1. -./data/pathway-genes.gmt.csv
2. -./data/unionpeak2.0
3. -peakscores from this [link](https://drive.google.com/file/d/10pjFdXjR_CqgZxTnUut7IceRpeDUXY2Y/view?usp=sharing)

**Step 2.** The output results will be saved inside folder '3'.

```shell
mkdir 3

R CMD BATCH '--args  starti=1 endi=295  MLmodel=3  outf="3/f-1-295.txt" TPRf="3/TPR-1-295.txt" FDRf="3/FDR-1-295.txt" predf="3/pred-1-295.txt" scoresf="3/scores-1-295.txt" toptf="3/tfs-1-295.txt"  cfmtx="3/cf-1-295.txt" toptf_score="3/tf_score-1-295.txt" toptf_cor="3/tf_cor-1-295.txt" '  functionML.R  &

```
**Step 3.** Run the  ./code/calculate_tpr_predict_genes.py script after *ensuring path to input files* (unionpeak2.0 and pathway-genes.gmt.csv) are correct.

```shell

python calculate_tpr_predict_genes.py

```
<!-- markdownlint-enable -->

Finally, output file names 'tpr_ml3_all_onto1.csv' will contain predicted genes for all the 295 gene-sets.


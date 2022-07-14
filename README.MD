## Introduction[![](./extras/images/pin.svg)](#introduction)

**GFPredict** is a R package. It predicts coding and non-coding genes related to (user input) gene-sets. GFPredict compares the similarity between the epigenomic and transcription factor binding patterns at the promoter sites of the input genes and other genes, to make predictions. [Click here](https://www.biorxiv.org/content/10.1101/2022.04.17.488570v1) to learn more about GFPredict.

---

## Table of contents[![](./extras/images/pin.svg)](#table-of-contents)
1. [GFPredict Package installation](#installation)

2. [GFPredict usage examples](#examples)

3. [Generate results used in the analyses of the paper](#results)
   - [Results of random forest algorithm for 295 ontological gene-sets](#ontology_results)
   - [PUBMED-based abstract mining](#pubmed_abtract_mining)
   - [Clustering of functions](#function_cluster)
   - [Validation of predicted lncRNAs in clusters using CRISPR screens](#cluster_crispr_screen)
   - [GFPredict testing using different CRISPR screens](#different_crispr_screen)


---

## 1. GFPredict Package installation[![](./extras/images/pin.svg)]<a name="installation"></a>

### Download zip file and install the package (recommended):

**Step 1.** Download zip file.

```shell
wget https://github.com/reggenlab/GFPredict/archive/refs/heads/main.zip

```
**Step 2.** Install the package in R session.

```r
devtools::install_local('main.zip')

```

### Or install through R commands.


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

### Generateing results of random forest algorithm for 295 ontological gene-sets in **3 steps:** <a name="ontology_results"></a>

**Step 1.** Clone the GitHub repository. 

```shell
#script
git clone https://github.com/reggenlab/GFPredict.git

```

**Step 2.** Run the below two commands and the output results will be saved inside directory '/3'. ***Make sure following R packages are installed:***

*  "ggplot2"
*  "randomForest" 
*  "glmnet"
*  "ROCR"
*  "xgboost"
*  "mltest"
*  "SHAPforxgboost"

**Note: It will take up to 30 minutes to completely run the code, check the progress opening functionalML.Rout file**

```shell
#move to the directory condtaining R script
cd GFPredict/code

#run the script functionML.R.
nohup R CMD BATCH '--args  starti=1 endi=295  MLmodel=3  outf="rf_result/f-1-295.txt" TPRf="rf_result/TPR-1-295.txt" FDRf="rf_result/FDR-1-295.txt" predf="rf_result/pred-1-295.txt" toptf="rf_result/tfs-1-295.txt"  cfmtx="rf_result/cf-1-295.txt" '  functionalML.R  &

```
**Step 3.** After getting output files in /3 directory, run the below command, predicted genes will be in the file.

```shell
#run this to get final predicted genes
python3 calculate_tpr_predict_genes.py

```

The output file 'predicted_genes.csv' will contain predicted genes for all the 295 gene-sets. Each row is a gene-set, corresponding columns contains predicted genes along with confidence score/tpr score


### PUBMED-based abstract mining <a name="pubmed_abtract_mining"></a>

For the predicted genes (file can also be found at ../data/predicted_genes) of random forest model: 
Run pubmed-based abstract mining and plot the number of occurrences of predicted gene terms and corresponding function terms.

```shell

python3 pubmed_abtract_mining.py

```


### Clustering of functions <a name="function_cluster"></a>

**Note:** R package 'Rtsne' and  'factoextra' must be installed to run this code.

Run the below code to get plot of cluster of functions based on shared top predictors.

```shell
Rscript clustering_functions.R
```


<!-- markdownlint-enable -->

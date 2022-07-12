
## Introduction[![](./extras/imgages/pin.svg)](#introduction)


# Installing GFPred R package

options(timeout=9999999)

devtools::install_github("reggenlab/GFPredict")

In case of difficulty in installation of libraries, "xgboost" and 'glmnet':

packageurl <- "http://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_0.90.0.2.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
install_version("glmnet", version = "3.0-2") #https://cran.r-project.org/src/contrib/Archive/glmnet/


# GFPredict 


Running the code to get predictions for different models.
Following R library packages are required:
1. ggplot2
2. randomForest 
3. caret  
4. e1071 
5. magrittr 
6. dplyr    
7. glmnet  
8. ROCR
9. xgboost 
10. mltest
11. SHAPforxgboost
12. tidyverse

Following python packages are required:
1. multiprocessing
2. concurrent.futures
3. itertools
4. glob

To get the predictions run the two scripts: functionML.R and calculate_tpr_predict_genes.py 
Running for 295 gene sets in parallel:

**Data:** Download the peak-count matrix (features) here: https://drive.google.com/drive/u/0/folders/1LSWhLrMwKxlrVfgIPR-ajCQnM4sBWozw 
Rest of the files are at the data directory 


**STEP 1:**
Running model 1 (Logistic Regression): 
  
  R CMD BATCH '--args  starti=1 endi=295  MLmodel=1  outf="1/f-1-295.txt" TPRf="1/TPR-1-295.txt" FDRf="1/FDR-1-295.txt" predf="1/pred-1-295.txt" scoresf="1/scores-1-295.txt" toptf="1/tfs-1-295.txt"  cfmtx="1/cf-1-295.txt"  '  functionML.R  &
  
  Running model 2 (Linear Regression): 
  
  R CMD BATCH '--args  starti=1 endi=295  MLmodel=2  outf="2/f-1-295.txt" TPRf="2/TPR-1-295.txt" FDRf="2/FDR-1-295.txt" predf="2/pred-1-295.txt" scoresf="2/scores-1-295.txt" toptf="2/tfs-1-295.txt"  cfmtx="2/cf-1-295.txt"  '  functionML.R  &
  
  Running model 3 (Random forest): 
  
  R CMD BATCH '--args  starti=1 endi=295  MLmodel=3  outf="3/f-1-295.txt" TPRf="3/TPR-1-295.txt" FDRf="3/FDR-1-295.txt" predf="3/pred-1-295.txt" scoresf="3/scores-1-295.txt" toptf="3/tfs-1-295.txt"  cfmtx="3/cf-1-295.txt" toptf_score="3/tf_score-1-295.txt" toptf_cor="3/tf_cor-1-295.txt" '  functionML.R  &
  
  Running model 4 (SVM):
  
  R CMD BATCH '--args  starti=1 endi=295  MLmodel=4  outf="4/f-1-295.txt" TPRf="4/TPR-1-295.txt" FDRf="4/FDR-1-295.txt" predf="4/pred-1-295.txt" scoresf="4/scores-1-295.txt" toptf="4/tfs-1-295.txt"  cfmtx="4/cf-1-295.txt"  '  functionML.R  &
  
  Running model 5 (XGBoost):
  
  R CMD BATCH '--args  starti=1 endi=295  MLmodel=4  outf="4/f-1-295.txt" TPRf="4/TPR-1-295.txt" FDRf="4/FDR-1-295.txt" predf="4/pred-1-295.txt" scoresf="4/scores-1-295.txt" toptf="4/tfs-1-295.txt"  cfmtx="4/cf-1-295.txt"  '  functionML.R  &


**STEP 2:**

Getting final predictions:
NOTE: The below script must be run 5 times by changing input files for 5 different model

nohup python calculate_tpr_predict_genes.py & (or python calculate_tpr_predict_genes.py)

Output: This will give output tpr_model3.to_csv('tpr_ml3_tfs_n_onto1.csv') which contains tpr values and predicted genes.





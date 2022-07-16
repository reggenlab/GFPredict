# This script takes in the assigned probablities by the random forest model
# It gives out final predicted genes for all the functions based on best confidence score value
import pandas as pd
import numpy as np

pathwaygenes = pd.read_csv('../data/pathway-genes.gmt.csv',header=None, sep="\t" ); #reads the ontology gene-sets files
unionpeaks = pd.read_csv('../data/unionpeak2.0',header=None, sep=" "); #curated list of promoters
pred_list = ['pred_scores1.csv', 'pred_scores2.csv', 'pred_scores3.csv', 'pred_scores4.csv', 
             'pred_scores5.csv', 'pred_scores6.csv', 'pred_scores7.csv', 'pred_scores8.csv',
             'pred_scores9.csv', 'pred_scores10.csv']
pred_i = []
for filename in pred_list:
    df = pd.read_csv(''.join(['..', '/data','/',filename]), index_col=None, header=0)
    pred_i.append(df)

pred_data = pd.concat(pred_i, axis=1, ignore_index=False)


def get_prbs_pgenes(row):          # Function to get the probabilities for each gene-set specified by index 'row'
    func = pathwaygenes.iloc[row,]
    pos_gene = func.dropna(axis=0)
    prbs = pred_data.iloc[:,row]
    prbs.index = unionpeaks.iloc[:,3]
    return pos_gene, prbs

def calculate_best_tpr(pos_gene, prbs):   # calculates best confidence score
    mat_tpr = pd.DataFrame()
    for t in np.arange(0, 1, 0.01): #Iterating between 0 and 1, to take threshold over probabilities
        t = round(t, 2)
        result = [index for (index, prbs) in enumerate(prbs) if prbs >= t]
        if len(result) != 0:
            pred_genes = unionpeaks.iloc[result, 3] #all the predicted genes above threshold t
            c_genes = pd.Series(list(set(pos_gene) & set(pred_genes)), dtype=pd.StringDtype())  # Take true positives from the predicted set
            tpr = len(c_genes) * 100 / len(set(pred_genes))  #Number of true positive genes divided by total number of predicted genes to get tpr/confidence score
            mat_tpr.loc[t, 1] = tpr
        max_tpr = max(mat_tpr.iloc[:,0])
        cutoff_max_tpr = mat_tpr.loc[mat_tpr.isin([max_tpr]).any(axis=1)].index[0]  # taking cutoff which yielded maximum tpr
        cut_genes = unionpeaks.iloc[[index for (index, prbs) in enumerate(prbs) if prbs >= cutoff_max_tpr], 3] # selecting genes with cutoff which yielded maximum tpr
        final_genes = pd.Series(list(set(cut_genes) - set(pos_gene)), dtype=pd.StringDtype()) #Taking final predicted genes after removing positives
    return max_tpr, cutoff_max_tpr, final_genes

def mk_mat(rangi):  
    mat_names = list() #creating list to name the columns of mat_names
    for i in range(1,2000):
        mat_names.append(str(i)+'_gene');
    mat_names[0:0] = ["func_name", "tpr", "cut_off"] #adding names to the list
    tpr_mat = pd.DataFrame(columns = mat_names); #empty dataframe to which calculated genes will be added
    for row in rangi:
        print("Number of gene-sets process:", row)
        pos_gene, prbs = get_prbs_pgenes(row) #taking positive genes and probabilities for the geneset indicated by index row
        if prbs.sum() != 0:
            max_tpr, cutoff_max_tpr, final_genes = calculate_best_tpr(pos_gene, prbs) #calculating the best tpr
            tpr_mat.loc[row,"func_name"] = pathwaygenes.iloc[row,0] #adding gene-set name to the dataframe
            tpr_mat.loc[row,"tpr"] = max_tpr #adding tpr to the dataframe
            tpr_mat.loc[row,"cut_off"] = cutoff_max_tpr #adding cutoff to the dataframe
            col_names = list()       #adding column names to the predicted genes (line 47 to 52)
            if len(final_genes) >= 2000:
                final_genes = final_genes[0:1999,]
            for j in range(1,len(final_genes) + 1):
                col_names.append(str(j)+'_gene')
            final_genes.index = col_names; 
            tpr_mat.loc[row, col_names] = final_genes # adding predicted genes to the dataframe 
    return tpr_mat

tpr_mat = mk_mat(range(0, 295)) #calling function mk_mat
tpr_mat.to_csv('predicted_genes.csv')

print("check 'predicted_genes.csv' file for the predicted genes")
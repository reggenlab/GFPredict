import pandas as pd
import csv
#from Bio import Entrez
import numpy as np
from os import listdir 
import glob
import os
import time
import random
from multiprocessing import Process, Pool
import concurrent.futures
from itertools import repeat
from networkx.algorithms.community.centrality import girvan_newman


pathwaygenes = pd.read_csv('~/omkar/non_coding/onto1/pathway-genes.gmt.csv',header=None, sep="\t" );
#pathwaygenes = pd.read_csv("~/omkar/non_coding/onto2/pathwaygenes2.0.csv", sep = "\t", header = None)

unionpeaks = pd.read_csv('~/omkar/non_coding/onto1/unionpeak2.0',header=None, sep=" ");
nc_genes = pd.read_csv('~/omkar/data/non_coding_rna_list.txt', header = 0, index_col = 0, sep = ",");
func_list1 = ['pred-1-295.txt', 'pred-296-590.txt', 'pred-591-885.txt', 'pred-886-1180.txt', 'pred-1181-1475.txt', 'pred-1476-1770.txt', 'pred-1771-2065.txt', 'pred-2066-2360.txt', 'pred-2361-2655.txt', 'pred-2656-2950.txt', 'pred-2951-3245.txt', 'pred-3246-3540.txt', 'pred-3541-3835.txt', 'pred-3836-4130.txt', 'pred-4131-4425.txt', 'pred-4426-4720.txt', 'pred-4721-5015.txt', 'pred-5016-5310.txt', 'pred-5311-5605.txt', 'pred-5606-5900.txt', 'pred-5901-5917.txt']
func_list2 = ['pred-1-295.txt', 'pred-296-590.txt', 'pred-591-885.txt', 'pred-886-1180.txt', 'pred-1181-1475.txt', 'pred-1476-1770.txt', 'pred-1771-2065.txt', 'pred-2066-2360.txt', 'pred-2361-2655.txt', 'pred-2656-2950.txt', 'pred-2951-3245.txt', 'pred-3246-3540.txt', 'pred-3541-3642.txt']

def get_prbs_pgenes(row):
    func = pathwaygenes.iloc[row,]
    pos_gene = func.dropna(axis=0)
    prbs = pred_data.iloc[:,row]
    prbs_shuff = prbs.sample(len(prbs))
    prbs.index = unionpeaks.iloc[:,3]
    prbs_shuff.index = unionpeaks.iloc[:,3]
    return pos_gene, prbs, prbs_shuff


def calculate_tpr_top(cut_genes, pos_gene, top, prbs):
    tpr_i = 0
    tpr_genes_unique = prbs.loc[cut_genes].groupby(prbs.loc[cut_genes].index).max()
    tpr_top_genes = tpr_genes_unique.sort_values(ascending=False)[0:top]
    c_genes = pd.Series(list(set(pos_gene) & set(tpr_top_genes.index)))
    if len(c_genes) != 0:
        tpr_i = len(c_genes) * 100 / len(tpr_top_genes)
    if top == 50:
        return tpr_i, tpr_top_genes.index
    return tpr_i

def calculate_best_tpr(pos_gene, prbs):
    mat_tpr = pd.DataFrame()
    for t in np.arange(0, 1, 0.01):
        t = round(t, 2)
        result = [index for (index, prbs) in enumerate(prbs) if prbs >= t]
        if len(result) != 0:
            pred_genes = unionpeaks.iloc[result, 3]
            c_genes = pd.Series(list(set(pos_gene) & set(pred_genes))) #Take true positives from the predicted set
            tpr = len(c_genes) * 100 / len(set(pred_genes))
            mat_tpr.loc[t, 1] = tpr
        max_tpr = max(mat_tpr.iloc[:,0])
        cutoff_max_tpr = mat_tpr.loc[mat_tpr.isin([max_tpr]).any(axis=1)].index[0]
        cut_genes = unionpeaks.iloc[[index for (index, prbs) in enumerate(prbs) if prbs >= cutoff_max_tpr], 3]
        final_genes = pd.Series(list(set(cut_genes) - set(pos_gene)))
        lnc_genes = np.intersect1d(final_genes, nc_genes)
        tpr_10 = calculate_tpr_top(cut_genes, pos_gene, 10, prbs)
        tpr_20 = calculate_tpr_top(cut_genes, pos_gene, 20, prbs)
        tpr_30 = calculate_tpr_top(cut_genes, pos_gene, 30, prbs)
        tpr_50, top50_genes = calculate_tpr_top(cut_genes, pos_gene, 50, prbs)
        tpr_tops = [tpr_10, tpr_20, tpr_30, tpr_50]
    return max_tpr, cutoff_max_tpr, final_genes, lnc_genes, tpr_tops, top50_genes 

def read_files(path, func_list):
    pred_file = pd.DataFrame()
    for i in func_list:
        words = [path, i]
        pathi = "/".join(words)
        print(pathi)
        file = pd.read_csv(pathi, sep = " ", header = 0)
        pred_file = pd.concat([pred_file, file], axis = 1)
    return pred_file

def mk_mat(rangi):
    mat_names = list()
    for i in range(1,2000):
        mat_names.append(str(i)+'_gene');
    mat_names[0:0] = ["func_name", "tpr", "cut_off", "tpr_10", "tpr_20", "tpr_30", "tpr_50"]
    tpr_mat = pd.DataFrame(columns = mat_names); tpr_mat_shuff = pd.DataFrame(columns = mat_names)
    lncr_mat = pd.DataFrame(columns = pathwaygenes.iloc[:,0], index = list(range(0,2500))); 
    lncr_mat_shuff = pd.DataFrame(columns = pathwaygenes.iloc[:,0], index = list(range(0,2500)))
    top50_mat = pd.DataFrame(columns = pathwaygenes.iloc[:,0], index = list(range(0,50))); 
    top50_mat_shuff = pd.DataFrame(columns = pathwaygenes.iloc[:,0], index = list(range(0,50)))
    for row in rangi:
        print(row)
        pos_gene, prbs, prbs_shuff = get_prbs_pgenes(row)
        if prbs.sum() != 0:
            max_tpr, cutoff_max_tpr, final_genes, lnc_genes, tpr_tops, top50_genes = calculate_best_tpr(pos_gene, prbs)
            max_tpr_shuff, cutoff_max_tpr_shuff, final_genes_shuff, lnc_genes_shuff, tpr_tops_shuff, top50_genes_shuff = calculate_best_tpr(pos_gene, prbs_shuff)
            tpr_mat.loc[row,"func_name"] = pathwaygenes.iloc[row,0]
            tpr_mat.loc[row,"tpr"] = max_tpr
            tpr_mat.loc[row,"cut_off"] = cutoff_max_tpr
            tpr_mat.loc[row,["tpr_10","tpr_20","tpr_30", "tpr_50"]] = tpr_tops
            tpr_mat_shuff.loc[row, "func_name"] = pathwaygenes.iloc[row,0]
            tpr_mat_shuff.loc[row,"tpr"] = max_tpr_shuff
            tpr_mat_shuff.loc[row,"cut_off"] = cutoff_max_tpr_shuff
            tpr_mat_shuff.loc[row,["tpr_10","tpr_20","tpr_30", "tpr_50"]] = tpr_tops_shuff
            top50_genes = pd.DataFrame(top50_genes); top50_genes.columns =[pathwaygenes.iloc[row,0]];
            top50_genes.index = list(range(0, len(top50_genes)))
            top50_genes_shuff = pd.DataFrame(top50_genes_shuff); top50_genes_shuff.columns = [pathwaygenes.iloc[row,0]];
            top50_genes_shuff.index = list(range(0, len(top50_genes_shuff)))
            top50_mat.loc[:,pathwaygenes.iloc[row,0]] = top50_genes
            top50_mat_shuff.loc[:,pathwaygenes.iloc[row,0]] = top50_genes_shuff
            col_names = list()
            if len(final_genes) >= 2000:
                final_genes = final_genes[0:1999,]
            for j in range(1,len(final_genes) + 1):
                col_names.append(str(j)+'_gene')
            final_genes.index = col_names; 
            col_names_shuff = list()
            if len(final_genes_shuff) >= 2000:
                final_genes_shuff = final_genes_shuff[0:1999,]
            for j in range(1,len(final_genes_shuff) + 1):
                col_names_shuff.append(str(j)+'_gene')
            final_genes_shuff.index = col_names_shuff
            lnc_genes = pd.DataFrame(lnc_genes); lnc_genes.columns = [pathwaygenes.iloc[row,0]]
            lnc_genes.index = list(range(0, len(lnc_genes)))
            lnc_genes_shuff = pd.DataFrame(lnc_genes_shuff); lnc_genes_shuff.columns = [pathwaygenes.iloc[row,0]]
            lnc_genes_shuff.index = list(range(0, len(lnc_genes_shuff)))
            tpr_mat.loc[row, col_names] = final_genes
            tpr_mat_shuff.loc[row, col_names_shuff] = final_genes_shuff
            lncr_mat.loc[:,pathwaygenes.iloc[row,0]] = lnc_genes
            lncr_mat_shuff.loc[:,pathwaygenes.iloc[row,0]] = lnc_genes_shuff
    return tpr_mat, lncr_mat, tpr_mat_shuff, lncr_mat_shuff, top50_mat, top50_mat_shuff


def call_calculation(filepath, func_list):
    global pred_data
    pred_data = read_files(filepath, func_list)    
    pred_data.columns = list(range(0, 5917))
    with concurrent.futures.ProcessPoolExecutor() as executor:
        rangi = [range(0, 600), range(601, 1200), range(1201, 1800), range(1801, 2400),
             range(2401,3000), range(3001, 3600), range(3601, 4200), range(4201, 4800), range(4801, 5400),
             range(5401, pathwaygenes.shape[0])]
        results_set = executor.map(mk_mat, rangi)
    sets = []
    for i in results_set:
        sets.append(i)
    tpr_set = pd.DataFrame(); lncr_set = pd.DataFrame(); tpr_set_shuff = pd.DataFrame();
    lncr_set_shuff = pd.DataFrame(); top50_set = pd.DataFrame(); top50_set_shuff = pd.DataFrame();
    for i in range(len(sets)):
        tpr_set = pd.concat([tpr_set, pd.DataFrame(sets[i][0])])
        lncr_set = pd.concat([lncr_set, pd.DataFrame(sets[i][1])])
        tpr_set_shuff = pd.concat([tpr_set_shuff, pd.DataFrame(sets[i][2])])
        lncr_set_shuff = pd.concat([lncr_set_shuff, pd.DataFrame(sets[i][3])])
        top50_set = pd.concat([top50_set, pd.DataFrame(sets[i][4])])
        top50_set_shuff = pd.concat([top50_set_shuff, pd.DataFrame(sets[i][5])])
    return tpr_set, lncr_set, tpr_set_shuff, lncr_set_shuff, top50_set, top50_set_shuff 


filepath3 = '~/omkar/data/non_coding_tfs_onto1/3'
tpr_model3, lncr_mat, tpr_model3_shuff, lncr_mat_shuff, top50_set, top50_set_shuff = call_calculation(filepath3, func_list1)
tpr_model3.to_csv('tpr_ml3_tfs_n_onto1.csv')
lncr_mat.to_csv('lnc_ml3_tfs_n_onto1.csv')
tpr_model3_shuff.to_csv('tpr_shuff_ml3_tfs_n_onto1.csv')
lncr_mat_shuff.to_csv('lnc_shuff_ml3_tfs_n_onto1.csv')
top50_set.to_csv('top50_ml3_tfs_n_onto1.csv')
top50_set_shuff.to_csv('top50_shuff_ml3_tfs_n_onto1.csv')
os.system("python3.7 ~/omkar/inform_me.py --text 'model3 onto1 tfs tpr is done'")



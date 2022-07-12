import pandas as pd
import csv
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


pathwaygenes = pd.read_csv('./data/pathway-genes.gmt.csv',header=None, sep="\t" ); 

unionpeaks = pd.read_csv('./data/unionpeak2.0',header=None, sep=" ");

func_list1 = ['pred-1-295.txt', 'pred-296-590.txt', 'pred-591-885.txt', 'pred-886-1180.txt', 'pred-1181-1475.txt', 'pred-1476-1770.txt', 'pred-1771-2065.txt', 'pred-2066-2360.txt', 'pred-2361-2655.txt', 'pred-2656-2950.txt', 'pred-2951-3245.txt', 'pred-3246-3540.txt', 'pred-3541-3835.txt', 'pred-3836-4130.txt', 'pred-4131-4425.txt', 'pred-4426-4720.txt', 'pred-4721-5015.txt', 'pred-5016-5310.txt', 'pred-5311-5605.txt', 'pred-5606-5900.txt', 'pred-5901-5917.txt']


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
    return max_tpr, cutoff_max_tpr, final_genes

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
    mat_names[0:0] = ["func_name", "tpr", "cut_off"]
    tpr_mat = pd.DataFrame(columns = mat_names);
    for row in rangi:
        print(row)
        pos_gene, prbs, prbs_shuff = get_prbs_pgenes(row)
        if prbs.sum() != 0:
            max_tpr, cutoff_max_tpr, final_genes = calculate_best_tpr(pos_gene, prbs)
            tpr_mat.loc[row,"func_name"] = pathwaygenes.iloc[row,0]
            tpr_mat.loc[row,"tpr"] = max_tpr
            tpr_mat.loc[row,"cut_off"] = cutoff_max_tpr
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
            tpr_mat.loc[row, col_names] = final_genes
    return tpr_mat


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
    tpr_set = pd.DataFrame();
    for i in range(len(sets)):
        tpr_set = pd.concat([tpr_set, pd.DataFrame(sets[i][0])])
    return tpr_set


filepath3 = '3/'
tpr_model3 = call_calculation(filepath3, func_list1)
tpr_model3.to_csv('tpr_ml3_tfs_n_onto1.csv')

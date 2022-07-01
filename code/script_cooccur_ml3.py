## This is the script to find the co-occurrence of predicted terms and input terms in PubMed abstracts

import pandas as pd
import glob
import csv
import numpy as np
from Bio import Entrez
from os import listdir 
import os
import re
from matplotlib import pyplot as plt
import random

pd.set_option('display.max_rows', 6000)
pd.set_option('display.max_columns', 200)
pd.set_option('display.width', 800)
pd.set_option('display.max_colwidth', None)

#pathwaygenes = pd.read_csv('~/omkar/non_coding/onto1/pathway-genes.gmt.csv',header=None, sep="\t" );
pathwaygenes = pd.read_csv("/home/omkarc/omkar/non_coding/onto2/pathwaygenes2.0.csv", sep = "\t", header = None, low_memory=False)
unionpeaks = pd.read_csv('/home/omkarc/omkar/non_coding/onto1/unionpeak2.0',header=None, sep=" ")


def text_process(func_word):     
    result = func_word[3:] #remove first three characters
    result = re.sub('[^A-Za-z0-9]+', ' ', result) #remove underscores
    result = re.sub(" \d+", "", result) # remove numbers 5536
    result = result.lower() 
    result = re.sub(r'(?:^| )\w(?:$| )', ' ', result).strip() #remove single letters
    stopwords = {'cell','of','small','in','is','he', 'to', 'from', 'by', 'on', 'the', 'or', 'like', 'layer',
                'ii', 'groups', 'into', 'type'}
    result  = ' '.join(filter(lambda x: x.lower() not in stopwords,  result.split()))
    if (len(result.split()) >= 2):
        stopwords1 = {'binding', 'protein', 'factor', 'activity', 'regulation', 'group', 'chemical', 'sensory', 'other', 'process', 
                    'species', 'positive', 'compound', 'cellular', 'particle', 'organism', 'involved', 'movement', 
                        'interaction', 'environment', 'pathway', 'signaling', 'coupled', 'mrna', 'response', 'negative',
                        'modified', 'response', 'left', 'right', 'formation', 'nucleotide', 'receptor', 'gene', 'complex',
                        'dependent', 'maintenance', 'process', 'acid'}
        result = ' '.join(filter(lambda x: x.lower() not in stopwords1,  result.split())) #second round
        words = result.split()
        result = " ".join(sorted(set(words), key=words.index)) #remove duplicate words
        result = re.sub(r'(?:^| )\w(?:$| )', ' ', result).strip() #remove single letters
    if (len(result.split()) >= 3):
        stopwords3 = {'dna'}
        result = ' '.join(filter(lambda x: x.lower() not in stopwords3,  result.split()))
    result = result.strip()
    intact_func = func_word
    result = ' OR '.join(result.split())
    return result, intact_func


func_names_list = list() # get list of the function names for usage
for funcs in pathwaygenes.iloc[:,0]:
    value, funci = text_process(funcs)
    func_names_list.append(funci)

#func_name to get the function name from main file & row to get the genes from data
def make_terms(data, func_name, row):
    func_word, intact_func = text_process(func_name) 
    gene_list = data.iloc[row,2:].dropna()
    gene_func_terms = list()
    for gene in gene_list:
        if (type(gene) == str):
            words = [func_word, "AND", gene]
            gene_func_terms.append(' '.join(words))
    return gene_func_terms, func_word, intact_func

def co_occurance(terms):
    Entrez.email = "omkarchandra.r@gmail.com"
    search_results = Entrez.read(
        Entrez.esearch(
            db="pubmed", term= terms, mindate=1990, maxdate=2021, datetype="pdat", usehistory="y", retmax='1000'
        )
    )
    count = int(search_results["Count"])
    pmids = list(search_results["IdList"])
    print(terms)
    print("counts:", count)
    print("pmids:", len(pmids))
    return count, pmids

def mk_mats():
    #master_count = pd.DataFrame(index = func_names_list)
    mat_names = list()
    for i in range(1,2000):
        mat_names.append(str(i)+'_counts');
    mat_names[0:0] = ["total_count", "func_name", "total_genes"]
    master_count = pd.DataFrame(index = func_names_list, columns = list(mat_names))
    #master_count.columns = [mat_names]
    pmid_mat = pd.DataFrame(columns = func_names_list, index = list(range(0,2000)))
    threshold_func = list()
    return master_count, pmid_mat, threshold_func

def call_func(data, thr):
    master_count, pmid_mat, threshold_func = mk_mats()
    for row in range(data.shape[0]):
        print("row_numer:", row)
        if data.iloc[row].loc["tpr"] >= thr and data.iloc[row].loc["tpr"] != 100: 
            func_count = 0
            func_array = []
            pmid_array = []
            func_name = data.iloc[row].loc["func_name"]
            gene_func_list, func_word, intact_func = make_terms(data, func_name, row)
            if len(gene_func_list) != 0:
                for terms in gene_func_list:
                    if len(str(gene_func_terms).split()) < 4:
                        continue
                    counts, pmids = co_occurance(terms)
                    func_count = +counts
                    func_array.append(counts)
                    pmid_array.append(pmids)
                if(sum(func_array) != 0):
                    threshold_func.append(intact_func)
                    master_count.loc[intact_func].loc["total_count"] = func_count
                    master_count.loc[intact_func].loc["func_name"] = func_name
                    master_count.loc[intact_func].loc["total_genes"] = len(data.iloc[row,4:].dropna())
                    func_array_col = []
                    for i in range(1, len(func_array) + 1):
                        func_array_col.append(str(i)+'_counts')
                    func_array = pd.DataFrame(func_array); func_array = func_array.transpose()
                    func_array.columns = [func_array_col]; func_array.index = [func_name]
                    master_count.loc[intact_func, func_array_col] = func_array.values
                    pmid_list = [item for sublist in pmid_array for item in sublist]
                    if(len(pmid_list) != 0):
                        pmid_list = pd.DataFrame(pmid_list); pmid_list.columns = [intact_func]
                        pmid_list.index = list(range(0, len(pmid_list)))
                        pmid_mat.loc[:, func_name] =  pmid_list
    master_count = master_count.loc[threshold_func] #selecting only those funcs above threshold
    return master_count, pmid_mat


def mak_null_freq(term_freq, thr):
    null_mat = pd.DataFrame(np.nan, index = term_freq.index, columns= range(term_freq.shape[1] + 10))
    for i in range(term_freq.shape[0]):
        print(i)
        sh_nm_func = term_freq.iloc[i].loc["func_name"]
        true_all = pathwaygenes.loc[pathwaygenes.iloc[:,0].str.contains(sh_nm_func),]
        true_genes = true_all.iloc[0,2:].dropna()
        all_genes = pd.Series(list(set(unionpeaks.iloc[:,3]) - set(true_genes)))
        nul_genes = random.sample(list(all_genes.iloc[:,]),  int(term_freq.iloc[i].loc['total_genes']))
        nul_genes = pd.DataFrame(nul_genes);  nul_genes = nul_genes.transpose();
        nul_genes.index = [sh_nm_func]
        null_mat.loc[sh_nm_func].iloc[0:len(nul_genes)] = ", ".join(nul_genes.iloc[0,:])
        null_mat.loc[sh_nm_func, "func_name"] = str(sh_nm_func)
        null_mat.loc[sh_nm_func, "tpr"] = thr
    print("ori mat:", term_freq.shape,"null mat:", null_mat.shape)
    null_freq, null_pmid = call_func(null_mat, thr)
    return null_freq, null_pmid

def plot_it(term_freq, null_freq):
    plot_mat = pd.DataFrame(index = term_freq.index, columns = ["c_p_f_pred", "c_p_f_null"])
    for i in range(term_freq.shape[0]):
        row_list = term_freq.iloc[i].dropna()
        counts_per_func = len(row_list[3:][ row_list[3:] != 0 ])
        plot_mat.loc[term_freq.index[i], "c_p_f_pred"] = int(counts_per_func)
        if i in range(null_freq.shape[0]):
            row_list_null = null_freq.iloc[i].dropna()
            counts_per_func_nul = len(row_list_null[3:])
            plot_mat.loc[null_freq.index[i], "c_p_f_null"] = int(counts_per_func_nul)
    plot_mat = plot_mat.fillna(0)
    plot_mat.to_csv("ml1plot_mat.csv")
    fig = plt.figure(figsize=(10,10))
    X = plot_mat.index
    predicted =  plot_mat.c_p_f_pred
    Random = plot_mat.c_p_f_null
    X_axis = np.arange(len(X))
    plt.bar(X_axis - 0.2, predicted, 0.4, label = 'Predicted')
    plt.bar(X_axis + 0.2, Random, 0.4, label = 'Random')
    plt.ylabel('Frequency of co-occurrence')
    plt.xlabel('gene sets')
    plt.title("Co-occurrence of ontology term with predicted gene terms, all-features, random forest, ontology 2")
    plt.legend()
    plt.savefig('ml3_all_onto2_tpr.png')

def final_call(data, thr):
    term_freq, pmid_mat = call_func(data, thr)
    null_freq, null_pmid = mak_null_freq(term_freq, thr)
    term_freq.to_csv("ml3term_freq_all_onto2_tpr.csv")
    pmid_mat.to_csv("ml3pmid_mat_all_onto2_tpr.csv")
    null_freq.to_csv("ml3null_freq_all_onto2_tpr.csv")
    null_pmid.to_csv("ml3null_pmid_all_onto2_tpr.csv")
    plot_it(term_freq, null_freq)

data = pd.read_csv('/home/omkarc/omkar/data/mk_barplots/tpr_plots/non_coding_all_n_onto2/3/lnc_ml3_all_onto2.csv', index_col = None, sep = ",")
final_call(data, thr = 60)


os.system("python3.7 ~/omkar/inform_me.py --text 'my selected pubmed is done'")


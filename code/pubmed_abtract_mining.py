import warnings
import pandas as pd
import numpy as np
from Bio import Entrez
import os
import re
from matplotlib import pyplot as plt
import random
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

pathwaygenes = pd.read_csv("../data/pathway-genes.gmt.csv", sep = "\t", header = None, low_memory=False)
unionpeaks = pd.read_csv('../data/unionpeak2.0',header=None, sep=" ")


def text_process(func_word):      # this function is to remove all the unwanted words
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

def co_occurance(terms): ## this function calls the pubmed api to go through abstracts to find the co-occurrence of predicted gene term and function term
    Entrez.email = "omkarchandra.r@gmail.com" #Feel free to change this to your peronal registered email address
    search_results = Entrez.read(
        Entrez.esearch(
            db="pubmed", term= terms, mindate=1990, maxdate=2021, datetype="pdat", usehistory="y", retmax='1000'
        )
    )
    count = int(search_results["Count"]) # The number of counts of the cooccurrence of gene term and function term
    pmids = list(search_results["IdList"]) # PMID of the articles found
    print(terms)
    print("counts:", count)
    print("pmids:", len(pmids))
    return count, pmids

def mk_mats(): #This function makes matrices in which the results will be saved
    mat_names = list()
    for i in range(1,2000):
        mat_names.append(str(i)+'_counts');
    mat_names[0:0] = ["number_evidence_found_each_gene", "func_name", "total_genes"]
    master_count = pd.DataFrame(index = func_names_list, columns = list(mat_names))
    pmid_mat = pd.DataFrame(columns = func_names_list, index = list(range(0,2000)))
    threshold_func = list()
    return master_count, pmid_mat, threshold_func

def call_func(data, thr):
    master_count, pmid_mat, threshold_func = mk_mats() #calling matrix making funcitons
    #for row in range(data.shape[0]):
    for row in range(36, data.shape[0]):
        print("row_numer:", row)
        if data.iloc[row].loc["tpr"] >= thr and data.iloc[row].loc["tpr"] != 100: 
            func_array = [] #to use the result to append to the main matrix
            pmid_array = [] #to use the result to 
            func_name = data.iloc[row].loc["func_name"]
            gene_func_list, func_word, intact_func = make_terms(data, func_name, row)
            print(gene_func_list)
            if len(gene_func_list) != 0:
                for terms in gene_func_list:
                    if (len(terms.split(' ')[0]) == 0): break
                    counts, pmids = co_occurance(terms)
                    func_array.append(counts)
                    pmid_array.append(pmids)
                if(sum(func_array) != 0):
                    threshold_func.append(intact_func)
                    master_count.loc[intact_func].loc["number_evidence_found_each_gene"] = sum(func_array)
                    master_count.loc[intact_func].loc["func_name"] = func_name
                    master_count.loc[intact_func].loc["total_genes"] = len(data.iloc[row,4:].dropna())
                    func_array_col = []
                    for i in range(1, len(func_array) + 1):
                        func_array_col.append(str(i)+'_counts')
                    func_array_df = pd.DataFrame(func_array); func_array_df = func_array_df.transpose()
                    func_array_df.columns = [func_array_col]; func_array_df.index = [func_name]
                    master_count.loc[intact_func, func_array_col] = func_array_df.values
                    pmid_list = [item for sublist in pmid_array for item in sublist]
                    if(len(pmid_list) != 0):
                        pmid_list = pd.DataFrame(pmid_list); pmid_list.columns = [intact_func]
                        pmid_list.index = list(range(0, len(pmid_list)))
                        pmid_mat[func_name] = pmid_list
    master_count = master_count.loc[threshold_func] #selecting only those funcs above threshold
    return master_count, pmid_mat


def mak_null_freq(term_freq, thr): #creating dummy data with random genes and running pubmed abstract for co-occurrence 
    null_mat = pd.DataFrame(np.nan, index=term_freq.index.values, columns=range(predicted_data.shape[1]) )
    null_mat.columns = predicted_data.columns
    for i in range(term_freq.shape[0]):
        print(i)
        sh_nm_func = term_freq.iloc[i].loc["func_name"]
        true_all = pathwaygenes.loc[pathwaygenes.iloc[:,0].str.contains(sh_nm_func),]
        true_genes = true_all.iloc[0,2:].dropna()
        all_genes = pd.Series(list(set(unionpeaks.iloc[:, 3]) - set(true_genes)), dtype=pd.StringDtype())
        nul_genes = random.sample(list(all_genes.iloc[:,]),  int(term_freq.iloc[i].loc['total_genes']))
        null_mat.loc[i, 'func_name'] = sh_nm_func
        null_array_col = []
        for ii in range(1, len(nul_genes) + 1):
             null_array_col.append(str(ii)+'_gene')
        null_mat.loc[i, null_array_col] = nul_genes
        null_mat.loc[i, "tpr"] = int(thr)
        null_mat.index = range(0, null_mat.shape[0])
    print("ori mat:", term_freq.shape,"null mat:", null_mat.shape)
    null_freq, null_pmid = call_func(null_mat, thr)
    return null_freq, null_pmid


def final_call(data, thr):
    term_freq, pmid_mat = call_func(data, thr)
    null_freq, null_pmid = mak_null_freq(term_freq, thr)
    term_freq.to_csv("term_freq.csv")
    pmid_mat.to_csv("pmid_mat.csv")
    null_freq.to_csv("null_freq.csv")
    null_pmid.to_csv("null_pmid.csv")

predicted_data = pd.read_csv('../data/predicted_genes.csv', index_col = None, sep = ",")
final_call(predicted_data, thr=60) #Calling all the functions and saving the data
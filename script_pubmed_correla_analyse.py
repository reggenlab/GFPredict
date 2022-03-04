### this is the script used for checking co-occurrence of the predicted terms and input terms in PubMed abstracts. 

import pandas as pd
import numpy as np
from Bio import Entrez
import os
import re

def text_process(func_word):     
    result = re.sub('[^A-Za-z0-9]+', ' ', func_word) #remove underscores
    result = re.sub(" \d+", "", result) # remove numbers 5536
    result = result.lower() 
    result = re.sub(r'(?:^| )\w(?:$| )', ' ', result).strip() #remove single letters
    stopwords = {'cell','of','small','in','is','he', 'to', 'from', 'by', 'on', 'the', 'or', 'like', 'layer',
                'ii', 'groups', 'into', 'type', 'containing', 'protein', 'receptor', 'organ'}
    result  = ' '.join(filter(lambda x: x.lower() not in stopwords,  result.split()))
    if (len(result.split()) >= 2):
        stopwords1 = {'binding', 'factor', 'activity', 'regulation', 'group', 'chemical', 'sensory', 'other', 'process', 
                    'species', 'positive', 'compound', 'cellular', 'particle', 'organism', 'involved', 'movement', 
                        'interaction', 'environment', 'pathway', 'signaling', 'coupled', 'mrna', 'response', 'negative',
                        'modified', 'response', 'left', 'right', 'formation', 'nucleotide', 'gene', 'complex',
                        'dependent', 'maintenance', 'process', 'acid'}
        result = ' '.join(filter(lambda x: x.lower() not in stopwords1,  result.split())) #second round
        words = result.split()
        result = " ".join(sorted(set(words), key=words.index)) #remove duplicate words
        result = re.sub(r'(?:^| )\w(?:$| )', ' ', result).strip() #remove single letters
    if (len(result.split()) >= 3):
        stopwords3 = {'dna', 'transporter', 'activity'}
        result = ' '.join(filter(lambda x: x.lower() not in stopwords3,  result.split()))
    result = result.strip()
    intact_func = func_word
    result = ' OR '.join(result.split())
    return result, intact_func


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

def make_terms(data, func_name, row):
    func_word, intact_func = text_process(func_name) 
    gene = data.iloc[row,0]
    gene_func_terms = list()
    words = [func_word, "AND", gene]
    gene_func_terms.append(' '.join(words))
    return gene_func_terms, func_word, intact_func

def call_func(data):
    count_list = list() 
    pmid_mat = pd.DataFrame()
    for row in range(data.shape[0]):
        func_name = data.iloc[row, 1] 
        gene_func_terms, func_word, intact_func = make_terms(data, func_name, row)
        if len(str(gene_func_terms).split()) < 4:
            continue 
        print(intact_func)
        counts, pmids = co_occurance(gene_func_terms)
        count_list.append(counts)
        pmid_mat = pd.concat([pmid_mat, pd.DataFrame(pmids)], axis = 1)
    return count_list, pmid_mat

correlation_analyse_pre = pd.read_csv("Correlation_Analyze.csv", sep = ",",  low_memory=False)
correlation_analyse = correlation_analyse_pre.dropna(axis=1, how='all')
all_gene = pd.concat([correlation_analyse.X1_gene, correlation_analyse.X2_gene, correlation_analyse.X3_gene
                      ,correlation_analyse.X4_gene])

all_funcs = pd.concat([correlation_analyse.X1_func, correlation_analyse.X2_func
                       , correlation_analyse.X3_func, correlation_analyse.X4_func])

data_pre = pd.concat([all_gene, all_funcs], axis = 1)
data = data_pre.dropna()
count_list, pmid_mat = call_func(data)

with open("count_result_correla_analyse.txt", "w") as output:
    output.write(str(count_list))
pmid_mat.to_csv("pmid_result_correla_analyse.csv")

os.system("python3.7 ~/omkar/inform_me.py --text 'correla_analyse is done'")


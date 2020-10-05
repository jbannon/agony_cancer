import networkx as nx
import pandas as pd
import numpy as np 
from scipy.stats import entropy
import math
import matplotlib.pyplot as plt
import fnmatch as fnm
import os
from collections import Counter
import sys

def sort_dict(dictionary):
	temp={}
	for i in sorted(dictionary.keys()):
		temp[i]=dictionary[i]
	return(temp)

def read_ranks(project, sample):
	base_path = "../data/results_data/tcga/"+project+"/"
	if sample =="normal":
		suffix = "_N"
	else:
		suffix = "_T"
	files = os.listdir(base_path)
	needs_unzip=True
	for f in files:
		if fnm.fnmatch(project+suffix+"_ranks.txt",f):
			needs_unzip=False
	if needs_unzip:
		print("needs unzip")
		uzip_string = "unzip "+base_path+project+suffix+"_ranks.zip -d "+base_path
		os.system(uzip_string)
	with open(base_path+project+suffix+"_ranks.txt","r") as istream:
		lines = istream.readlines()
	rank_dict = {} # keys are ranks, values is a list of ids of genes
	rank_multiset = []
	for line in lines:
		temp = line.rstrip().split()
		rank = int(temp[1])
		gene = int(temp[0])
		rank_multiset.append(rank)
		if rank in rank_dict.keys():
			gene_list = rank_dict[rank]
			gene_list = gene_list +[gene]
			rank_dict[rank]= gene_list
		else:
			rank_dict[rank] = [gene]
	rank_dict = sort_dict(rank_dict)
	return rank_dict, rank_multiset

def get_gene_names(d,gn):
	n ={}
	for k in d.keys():
		ids = d[k]
		n[k]=gn.loc[ids,'external_gene_name'].tolist()
	return n 

def count_ts_at_ranks(d,ts_genes):
	names = {}
	counts = {}
	for k in d.keys():
		genes = d[k]
		dx = [g in ts_genes for g in genes]
		gdx = [i for i in range(len(dx)) if dx[i]]
		counts[k]=len(gdx)
		names[k]=[genes[i] for i in gdx]
	return names, counts



def rank_analysis(project, normal_ranks, tumor_ranks,ts_genes):
	## needs to: 
	
	
	

	## 	 compute the distinct number of ranks needed for tumor and normal
	nrc = len(normal_ranks.keys())
	trc = len(tumor_ranks.keys())

	##   compute number of genes at each rank
	gene_map = pd.read_csv("../data/input_data/tcga/"+project+"/index_gene_map.csv")
	
	##   compute the number of TS genes at each rank
	normal_named = get_gene_names(normal_ranks,gene_map)
	tumor_named = get_gene_names(tumor_ranks, gene_map)

	norm_name,norm_count = count_ts_at_ranks(normal_named,ts_genes)
	tumr_name, tumr_count = count_ts_at_ranks(tumor_named,ts_genes)

	nranks = max(max(normal_ranks.keys()),max(tumor_ranks.keys()))
	
	df = pd.DataFrame({'ranks':range(nranks+1), 'tumor':[0 for i in range(nranks+1)],'normal':[0 for i in range(nranks+1	)]})
	for k in norm_count.keys():
		df.iloc[k,2]=norm_count[k]
	for k in tumr_count.keys():
		df.iloc[k,1]=tumr_count[k]
	df.to_csv("../data/results_data/tcga/"+project+"/ts_counts.csv")
	
	df = pd.DataFrame({'ranks':range(nranks+1), 'tumor':[0 for i in range(nranks+1)],'normal':[0 for i in range(nranks+1	)]})
	for k in normal_ranks.keys():
		df.iloc[k,2]=len(normal_ranks[k])
	for k in tumor_ranks.keys():
		df.iloc[k,1]=len(tumor_ranks[k])

	df.to_csv("../data/results_data/tcga/"+project+"/rank_counts.csv")
	


if __name__ == '__main__':
	with open("../data/tumor_suppressors.txt") as istream:
		ts=istream.readlines()
	ts = [x.rstrip() for x in ts]
	nrc, trc = 0,0
	nrcd, trcd = {},{}
	projects=["BLCA","BRCA", "COAD", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC","PRAD", "STAD", "THCA"]
	for project in projects:
		normal_ranks,normal_multiset = read_ranks(project, "normal")
		tumor_ranks,tumor_multiset = read_ranks(project, "tumor")
		rank_analysis(project,normal_ranks, tumor_ranks, ts)


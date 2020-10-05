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
	return rank_dict, rank_multiset

def compute_entropy(ranklist):
	value,counts = np.unique(ranklist, return_counts=True)
	return(entropy(counts/np.sum(counts)))

def plot_rank_distribution(project, sample, rank_list,log=False,entropy=True):
	ranks= Counter(rank_list).keys()
	counts = Counter(rank_list).values()
	if log:
		counts = [math.log(x) for x in list(counts)]
		title_string = project + " "+sample +" Rank Barchart (log)"
		fname = "../figs/tcga/"+project+"/"+sample+"_rank_distro_log.png"
		ylab = "log count"
	else:
		title_string = project + " "+sample +" Rank Barchart"
		fname = "../figs/tcga/"+project+"/"+sample+"_rank_distro.png"
		ylab = "Count"
	ax = plt.subplot(111)
	plt.bar(ranks, counts)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	plt.title(title_string)
	plt.xlabel("Rank")
	plt.xticks(np.arange(len(ranks)),ranks)
	plt.ylabel(ylab)
	plt.savefig(fname)
	plt.close()
	if entropy:
		H = compute_entropy(rank_list)
		with open("../figs/tcga/"+project+"/"+project+"_"+sample+"entropy.txt","w") as of:
			of.writelines(str(H))

def fetch_global_agonies(projects):
	normal_agonies = []
	tumor_agonies = []
	for project in projects:
		with open("../data/results_data/tcga/"+project+"/"+project+"_N_total.txt") as f:
			lines = f.readlines()
		normal_agonies.append(float(lines[0]))
		with open("../data/results_data/tcga/"+project+"/"+project+"_T_total.txt") as f:
			lines = f.readlines()
		tumor_agonies.append(float(lines[0]))
	return pd.DataFrame({"cancer":projects, "normal_agony":normal_agonies, "log_normal_agonies":[math.log(z) for z in normal_agonies],
	 "tumor_agonies":tumor_agonies, "log_tumor_agonies":[math.log(z) for z in tumor_agonies]})

def differential_rank_analysis(project, normal_ranks, tumor_ranks,k=0):
	normal_genes = set(normal_ranks[k])
	tumor_genes = set(tumor_ranks[k])
	both_samples = list(normal_genes.intersection(tumor_genes))
	tumor_genes = list(tumor_genes.difference(normal_genes))

	gene_map = pd.read_csv("../data/input_data/tcga/"+project+"/index_gene_map.csv")
	tumor_high_rank = gene_map.loc[tumor_genes,'external_gene_name'].tolist()
	tumor_high_rank=[x+"\n" for x in tumor_high_rank]
	with open("../data/results_data/tcga/"+project+"/tumor_only_at_rank_"+str(k)+".txt","w") as of:
		of.writelines(tumor_high_rank)

	normal_genes = list(normal_genes.difference(tumor_genes))
	normal_high_rank = gene_map.loc[normal_genes,'external_gene_name'].tolist()
	normal_high_rank=[x+"\n" for x in normal_high_rank]
	with open("../data/results_data/tcga/"+project+"/normal_only_at_rank_"+str(k)+".txt","w") as of:
		of.writelines(normal_high_rank)

	both_genes = gene_map.loc[both_samples,'external_gene_name'].tolist()
	both_genes=[x+"\n" for x in both_genes]
	with open("../data/results_data/tcga/"+project+"/both_at_rank_"+str(k)+".txt","w") as of:
		of.writelines(both_genes)


def count_at_ranks(rdict):
	counts = {}
	for key in rdict.keys():
		counts[key]=len(rdict[key])
	return counts

def make_rank_table(nr,tr):
	ncounts = count_at_ranks(nr)
	tcounts = count_at_ranks(tr)
	max_rank = max(max(ncounts.keys()),tcounts.keys())



def plot_network(graph):
	G = nx.read_weighted_edgelist("../data/results_data/tcga/LIHC/LIHC_N_edges.txt")
	pos = nx.circular_layout(G)
	f = plt.figure()
	nx.draw(G, pos=pos, ax=f.add_subplot(111))
	f.savefig("./graph.png")

def read_edgelist(project):
	G = nx.read_weighted_edgelist(path)
	return(G)

if __name__ == '__main__':
	projects=["BLCA","BRCA", "COAD", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC","PRAD", "STAD", "THCA"]
	global_ag = fetch_global_agonies(projects)
	global_ag.to_csv("../data/results_data/tcga/ALL/global_agonies.csv")	
	for project in projects:
		normal_ranks,normal_multiset = read_ranks(project, "normal")
		tumor_ranks,tumor_multiset = read_ranks(project, "tumor")
		differential_rank_analysis(project,normal_ranks, tumor_ranks)
		norm_counts = count_at_ranks(normal_ranks)
		tumor_counts = count_at_ranks(tumor_ranks)
		print(project)
		print(norm_counts)
		print(tumor_counts)
		#plot_rank_distribution(project,"Normal",normal_multiset)
		#plot_rank_distribution(project,"Normal",normal_multiset,log=True,entropy=False)
		#plot_rank_distribution(project,"Tumor",tumor_multiset)
		#plot_rank_distribution(project,"Tumor",tumor_multiset,log=True,entropy=False)

	
import os
import sys

def fetch_counts(path):
	
	with open(path+"tdims.txt",'r') as f:
		lines = f.readlines()
	
	tsamples = lines[0]
	tgenes = lines[1]
	
	with open(path+"ndims.txt",'r') as f:
		lines = f.readlines()
	nsamples = lines[0]
	ngenes = lines[1]
	return [int(tsamples), int(nsamples)], [int(tgenes), int(ngenes)]


def replace_rename(ifname, ofname):
	with open(ifname,"r") as f:
		lines = f.readlines()
	for i in range(len(lines)):
		lines[i]=lines[i].replace(","," ")
	with open(ofname,"w") as of:
		of.writelines(lines)

def edit_fixer_file(cancer,n_genes, n_samples,sample_type, n_boot=10):
	if sample_type=="normal":
		s="_n"
	else:
		s="_t"
	with open('pphi_bs.c','r') as f:
		lines = f.readlines()

	lines[8] = '#define NROW '+str(n_genes)+'\n'
	lines[9] = '#define TSAMPLE_COUNT '+str(n_samples)+'\n'
	lines[10] = '#define NBOOTSTRAP '+str(n_boot)+'\n'

	with open('phixer_'+str(cancer)+s+'.c','w') as of:
		of.writelines(lines)
	os.system("gcc -Wall phixer_"+cancer+s+".c"+" -fopenmp -o phixer_"+cancer+s+".out")
	return "phixer_"+cancer+s+".out"

def main(cancer,sample):
	
	path = "../data/input_data/tcga/"+cancer+"/"
	sample_counts, gene_counts = fetch_counts(path)

	if sample == "normal":
		idx=1
		suffix = "_N"
	else:
		idx=0
		suffix = "_T"
	
	phixer_call = edit_fixer_file(cancer,n_genes = gene_counts[idx], n_samples=sample_counts[idx],sample_type=sample)
	phixer_call = phixer_call+" ../data/input_data/tcga/"+cancer+"/"+sample+"_expression.txt"
	os.system(phixer_call)
	phixer_out_string = "pruned_*_"+cancer+"_"+sample+"_expression.txt"
	new_edgelist_string = cancer+suffix+"_edges.txt"
	replace_rename(phixer_out_string,new_list_string)
	agony_ostring = cancer+suffix+"_ranks.txt"
	os.system("agony_weighted/agony -i "+new_edgelist_string+" -o "+agony_ostring+" -w")
	

if __name__=='__main__':
	iname = sys.argv[1]
	oname =sys.argv[2]
	main(cancer, sample)
	

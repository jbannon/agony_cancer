import sys
import os
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

def edit_fixer_file(cancer,n_genes, n_samples, n_boot=10):
	with open('pphi_bs.c','r') as f:
		lines = f.readlines()

	lines[8] = '#define NROW '+str(n_genes)+'\n'
	lines[9] = '#define TSAMPLE_COUNT '+str(n_samples)+'\n'
	lines[10] = '#define NBOOTSTRAP '+str(n_boot)+'\n'

	with open('phixer_'+str(cancer)+'.c','w') as of:
		of.writelines(lines)
	os.system("gcc -Wall phixer_"+cancer+".c"+" -fopenmp -o phixer_"+cancer+".out")


def main(cancer,sample):
	path = "../data/input_data/tcga/"+cancer+"/"
	sample_counts, gene_counts = fetch_counts(path)
	if sample == "normal":
		idx=1
	else:
		idx=0
	edit_fixer_file(cancer,n_genes = gene_counts[idx], n_samples=sample_counts[idx])
	


if __name__ == '__main__':
	cancer = sys.argv[1]
	sample = sys.argv[2]
	main(cancer, sample)

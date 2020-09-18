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
	return [tsamples, nsamples], [tgenes, ngenes]

def edit_fixer_file(cancer,sample_type, n_genes, n_samples, n_boot=10):
	with open('pphi_bs.c','r') as f:
		lines = f.readlines()

	lines[8] = '#define NROW '+str(n_genes)+'\n'
	lines[9] = '#define TSAMPLE_COUNT '+str(n_samples)+'\n'
	lines[10] = '#define NBOOTSTRAP '+str(n_boot)+'\n'

	with open('phixer_'+str(cancer)+'_'+sample_type+'.c','w') as of:
		of.writelines(lines)
	return('phixer_'+str(cancer)+'_'+sample_type+'.c')

def run_phixer(cancer,sample_type, n_genes, n_samples, n_boot=10):
	print(n_genes)
	print(n_samples)
	print(n_boot)
	expr_file = "../data/input_data/tcga/"+cancer+"/"
	if sample_type=="t":
		expr_file = expr_file+"tumor_expression.txt"
	else:
		expr_file = expr_file +"normal_expression.txt"

	phixer_file = edit_fixer_file(cancer,sample_type, n_genes, n_samples, n_boot=10)
	phixer_string = "phixer_"+cancer+"_"+sample_type+".out"
	os.system("gcc -Wall "+phixer_file+" -fopenmp -o "+phixer_string) # compile phixer for this data
	os.system("ulimit -s 1300000000")
	os.system("export GOMP_STACKSIZE=2000000")
	os.system("./"+phixer_string+" "+expr_file) # run phixer on the file
	
	# phixer outputs into the same directory as the .out file 
	# so shoud have a follow-up command to make move

def main(cancer):
	path = "../data/input_data/tcga/"+cancer+"/"
	sample_counts, gene_counts = fetch_counts(path)
	for sample_type,n_samples,n_genes in zip(["t","n"],sample_counts, gene_counts):
		run_phixer(cancer, sample_type=sample_type, n_genes=n_genes, n_samples=n_samples)
	


if __name__ == '__main__':
	cancer = sys.argv[1]
	main(cancer)



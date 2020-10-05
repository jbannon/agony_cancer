import os
import sys
import fnmatch

def threshold_edges(file, cutoff):
	print("not implemented yet")
	""" matlab thresholding function from phixer

	function thresholded_nw = phixer_threshold(pruned_file, threshold)
    if nargin < 1
        error('Missing input file\n');
    end
    % Default threshold -- keep edges whoes weight is in top 1 percentile.
    if nargin < 2
        threshold = 0.99;
    end
    
    % Load input file
    thresholded_nw =  dlmread(pruned_file, ',');
    
    % Convert to matrix
    thresholded_nw = spconvert(thresholded_nw);
    
    % Drop threshold percentile edges and compute total degree
    x = nonzeros(thresholded_nw); 
    t = quantile(x, threshold); 
    thresholded_nw(thresholded_nw<t) = 0;
    
    % Print number of strongly connected components
    fprintf('Graph SCC: %d\n',  graphconncomp(thresholded_nw));

end"""

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
	files=os.listdir(".")
	temp=""
	for f in files:
		if fnmatch.fnmatch(f,ifname):
			temp=f
	ifname=temp
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
	print("gcc -Wall phixer_"+cancer+s+".c"+" -fopenmp -o phixer_"+cancer+s+".out")
	return "./phixer_"+cancer+s+".out"

def main(cancer,sample,zip_files=True, threshold=False, cutoff=0.99):
	
	path = "../data/input_data/tcga/"+cancer+"/"
	sample_counts, gene_counts = fetch_counts(path)

	if sample == "normal":
		idx=1
		suffix = "_N"
	else:
		idx=0
		suffix = "_T"
	
	phixer_call = edit_fixer_file(cancer,n_genes = gene_counts[idx], n_samples=sample_counts[idx],sample_type=sample)
	phixer_call = phixer_call+" ../data/input_data/tcga/"+cancer+"/"+cancer+"_"+sample+"_expression.txt"
	os.system("ulimit -s 1300000000")
	os.system("export GOMP_STACKSIZE=2000000")
		
	os.system(phixer_call)
	phixer_out_string = "pruned_*_"+cancer+"_"+sample+"_expression.txt"
	new_edgelist_string = cancer+suffix+"_edges.txt"
	replace_rename(phixer_out_string,new_edgelist_string)
	agony_ostring = cancer+suffix+"_ranks.txt"
	os.system("agony_weighted/agony -i "+new_edgelist_string+" -o "+agony_ostring+" -w")
	os.system("mv "+"phixer_"+cancer+suffix.lower()+".c"+" ./code_cache/"+cancer+"_code/")
	os.system("mv "+"phixer_"+cancer+suffix.lower()+".out"+" ./code_cache/"+cancer+"_code/")
	os.system("mv "+cancer+suffix+".sh"+" ./code_cache/"+cancer+"_code/")
	if zip_files:
		os.system("zip "+cancer+suffix+"_ranks.zip "+cancer+suffix+"_ranks.txt ")
		os.system("zip "+cancer+suffix+"_edges.zip "+cancer+suffix+"_edges.txt ")
		os.system("mv "+cancer+suffix+"_ranks.zip " +"../data/results_data/tcga/"+cancer+"/")
		os.system("mv "+cancer+suffix+"_edges.zip " +"../data/results_data/tcga/"+cancer+"/")	
		os.system("mv "+cancer+suffix+"_total.txt " +"../data/results_data/tcga/"+cancer+"/")
		os.system("rm "+cancer+suffix+"_ranks.txt")
		os.system("rm "+cancer+suffix+"_edges.txt")
	else:
		print("not yet implemented")

if  __name__=='__main__':
	cancer = sys.argv[1]
	sample =sys.argv[2]
	print(cancer)
	print(sample)
	sys.exit(0)
	main(cancer, sample)
	

import os
import sys

def replace_rename(ifname, ofname):
	with open(ifname,"r") as f:
		lines = f.readlines()
	for i in range(len(lines)):
		lines[i]=lines[i].replace(","," ")
	with open(ofname,"w") as of:
		of.writelines(lines)


def main():
	os.system("gcc -Wall pphi_bs.c -fopenmp -o phixer_.out")
	os.system("./phixer.out test_100_585.txt")
	replace_rename("pruned_nw_bs10_bins12_test_100_585.txt","pruned_nw_bs10_bins12_test_100_585.txt")

	os.system("agony_weighted/agony -i pruned_nw_bs10_bins12_test_100_585.txt -o res.txt -w")

if __name__=='__main__':
	iname = sys.argv[1]
	oname =sys.argv[2]
	main()
	

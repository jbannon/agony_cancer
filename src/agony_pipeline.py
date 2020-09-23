import os



os.system("gcc -Wall pphi_bs.c -fopenmp -o phixer_.out")
os.system("./phixer.out test_100_585.txt")
os.system("agony_weighted/agony -i pruned*.txt")
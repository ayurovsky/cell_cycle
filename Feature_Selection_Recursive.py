import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json
import sys

prefix = "SC_df_z_train"
#prefix = "microarray_leg_test"

canonical_gene_list = list() 
with open(prefix + '_genes.csv', 'r') as f:
	reader = csv.reader(f)
	canonical_gene_list = list(reader)[0]
canonical_gene_list_indeces = list(range(0, len(canonical_gene_list)))

gene_search_list = canonical_gene_list_indeces.copy()
gene_idx_best_so_far_list = list()
gene_name_best_so_far_list = list()
accuracies_list = list()
print(len(gene_search_list))
k = 10
while(len(gene_idx_best_so_far_list) < len(canonical_gene_list)): #52):
	counter = 0
	new_genes_dict = dict()

	max_jobs = 500
	for i in range(0, int(len(gene_search_list)/max_jobs)+1):
		processes = []
		for m in range(0, max_jobs):
			j = i*max_jobs + m
			if (j == len(gene_search_list)):
				break
			args = ['./Feature_Selection_Recursive', str(k), prefix]
			args += gene_idx_best_so_far_list
			args.append(str(gene_search_list[j]))
			p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)
			processes.append(p)

		for p in processes:
			stdout, stderr = p.communicate()
			#print(stderr.decode("utf-8"))
			#print(canonical_gene_list[gene_search_list[counter]])
			val = stdout.decode("utf-8")
			#print(val)
			new_genes_dict[gene_search_list[counter]] = val 
			counter += 1

	ListOfSortedGeneTuples = sorted(new_genes_dict.items(), reverse=True, key=lambda v : v[1])
	accuracies_list.append(int(ListOfSortedGeneTuples[0][1]))
	gene_idx_best_so_far_list.append(str(ListOfSortedGeneTuples[0][0]))
	gene_name_best_so_far_list.append(canonical_gene_list[ListOfSortedGeneTuples[0][0]])
	gene_search_list.remove(ListOfSortedGeneTuples[0][0])

	print("Iteration over")
	print(gene_name_best_so_far_list)
	print(accuracies_list)
	sys.stdout.flush()

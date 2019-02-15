import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json


canonical_gene_list = list() 
canonical_gene_list_indeces = list(range(0, 252))
with open('SC_df_z_train_genes.csv', 'r') as f:
	reader = csv.reader(f)
	canonical_gene_list = list(reader)[0]

gene_search_list = canonical_gene_list_indeces.copy()
gene_idx_best_so_far_list = list()
gene_name_best_so_far_list = list()
accuracies_list = list()
k = 7
while(len(gene_idx_best_so_far_list) < 252): #52):
	processes = []
	for j in range(0, len(gene_search_list)): 
		args = ['./Feature_Selection_Recursive', str(k)]
		args += gene_idx_best_so_far_list
		args.append(str(gene_search_list[j]))
		p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)
		processes.append(p)

	counter = 0
	new_genes_dict = dict()
	for p in processes:
		stdout, stderr = p.communicate()
		#print(stderr.decode("utf-8"))
		#print(canonical_gene_list[gene_search_list[counter]])
		val = stdout.decode("utf-8")
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

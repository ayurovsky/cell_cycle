import pandas as pd
import csv
import sys
import scipy.spatial.distance
from scipy.spatial.distance import euclidean
import json

SC_df = pd.DataFrame()

column_names = []
with open('Dataset930.csv') as f:
	reader = csv.reader(f, delimiter="\t")
	next(reader)
	counter = 0
	full_row = []
	gene_name = ""
	for row in reader:
		if (counter < 930):
			phase = ""
			if (row[1] == "g0/g1"):
				phase = "G1"
			elif (row[1] == "s"):
				phase = "S"
			elif (row[1] == "g2/m"):
				phase = "G2M"
			name = "cell_" + str(counter) + "_" + phase
			column_names.append(name)
		if (counter%930 == 0 and counter > 0):
			SC_df[gene_name] = full_row
			full_row = []
		full_row.append(2**float(row[9]))
		gene_name = row[6]
		counter += 1
SC_df = SC_df.T
SC_df.columns = column_names
print(SC_df.index)
print(SC_df.shape)
sys.stdout.flush()


# split into test and train
SC_df_train = SC_df.iloc[:,:500]
print(SC_df_train.shape)
num_s = 0
num_g2 = 0
num_g1 = 0
for cell_name in SC_df_train:
	if ("_S" in cell_name):
		num_s += 1
	elif ("_G2" in cell_name):   
		num_g2 += 1
	elif ("_G1" in cell_name):
		num_g1 += 1
print("num S: " + str(num_s))
print("num G2: " + str(num_g2))
print("num G1: " + str(num_g1))
SC_df_test = SC_df.iloc[:,500:]
print(SC_df_test.shape)
num_s = 0
num_g2 = 0
num_g1 = 0
for cell_name in SC_df_test:
	if ("_S" in cell_name):
		num_s += 1
	elif ("_G2" in cell_name):   
		num_g2 += 1
	elif ("_G1" in cell_name):
		num_g1 += 1
print("num S: " + str(num_s))
print("num G2: " + str(num_g2))
print("num G1: " + str(num_g1))
sys.stdout.flush()

# transform raw data to log1p before taking z-scores
import numpy as np
from scipy import stats

SC_df_train_log1p = SC_df_train.apply(np.log1p)
z_arr = stats.zscore(SC_df_train_log1p, axis=1)
SC_df_z_train = pd.DataFrame(z_arr, index = SC_df_train.index, columns = SC_df_train.columns)

SC_df_test_log1p = SC_df_test.apply(np.log1p)
z_arr = stats.zscore(SC_df_test_log1p, axis=1)
SC_df_z_test = pd.DataFrame(z_arr, index = SC_df_test.index, columns = SC_df_test.columns)

# # v1.3 k-NN weighted by 1/d
# do simple nearest neighbor classifier for k = 1...10
# for now try all genes
import scipy.spatial.distance
from scipy.spatial.distance import euclidean


def run_process(solved_gene_list, gene_list, run_index):
	#gene_list = ['ALAS1', 'ABCF1']
	#run_index = 0

	total_tested_dict = dict()
	total_correct_dict = dict()
	for gene_name in gene_list:
		total_tested_dict[gene_name] = 0
		total_correct_dict[gene_name] = 0

	#REMOVE RANGE - doing accuracy for only 20 TEST cells
	#cell_test_counter = 0
	for cell_test_name in SC_df_z_test: # go through all the test cells, compare each seperately to the "train" cells
		#cell_test_counter += 1
		#if (cell_test_counter > 20):
		#	break
		# get the list of distances for each of the cells
		dist_values_dict = dict()
		for gene_name in gene_list:
			dist_values_dict[gene_name] = list()
		#print(cell_test_name)
		sys.stdout.flush()
		for cell_train_name in SC_df_z_train:
			for gene_name in gene_list:
				proposed_gene_list = solved_gene_list.copy() + [gene_name]
				v_1 = SC_df_z_test.loc[proposed_gene_list,cell_test_name].values
				v_2 = SC_df_z_train.loc[proposed_gene_list,cell_train_name].values
				dist = euclidean(v_1, v_2)
				dist_values_dict[gene_name].append((cell_train_name, dist))
		# now sort them based on distance
		for gene_name in gene_list: 
			dist_values_dict[gene_name].sort(key=lambda x: x[1])
			# now get accuracies for k-nearest neighbor, k = 5
			m = 5
			total_G2M = 0
			total_G1 = 0
			total_S = 0
			for k in range(0,m):
				k_cell_train_name = dist_values_dict[gene_name][k][0]
				distance = dist_values_dict[gene_name][k][1]
				if ("_G2M" in k_cell_train_name):
					total_G2M += 1/distance
				elif ("_G1" in k_cell_train_name):
					total_G1 += 1/distance
				elif ("_S" in k_cell_train_name):
					total_S += 1/distance
			if (total_G2M == max([total_G2M, total_G1, total_S])):
				k_cell_train_name = "train_G2M"
			elif (total_G1 == max([total_G2M, total_G1, total_S])):
				k_cell_train_name = "train_G1"
			elif (total_S == max([total_G2M, total_G1, total_S])):
				k_cell_train_name = "train_S"
			if (("_G2M" in cell_test_name and "_G2M" in k_cell_train_name) or ("_G1" in cell_test_name and "_G1" in k_cell_train_name) or ("_S" in cell_test_name and "_S" in k_cell_train_name)):
				total_correct_dict[gene_name] += 1
			total_tested_dict[gene_name] += 1

	#for gene_name in gene_list:
	#	print("Accuracy for " + gene_name + " is " + str(total_correct_dict[gene_name]/total_tested_dict[gene_name]))

	json_name = "gene_list_len_" + str(len(solved_gene_list)) + "_run_index_" + str(run_index) + ".json"
	with open(json_name, "w") as write_file:
		json.dump(total_correct_dict, write_file)

from multiprocessing import Process
import sys
import os
import re

gene_search_list = SC_df_z_train.index.tolist()
gene_best_so_far_list = list()
gene_best_num_correct = 0 

while(len(gene_best_so_far_list) < 250):
	processes = []
	max_processors = 60
	num_per_process = max(1,int(len(gene_search_list)/max_processors))
	process_counter = 1 
	for start in range(0,len(gene_search_list),num_per_process):
		print("starting process " + str(process_counter))
		process_counter += 1
		stop = min(len(gene_search_list),start + num_per_process)
		p = Process(target=run_process, args=(gene_best_so_far_list, gene_search_list[start:stop], start))
		p.start()
		processes.append(p)

	for p in processes:
		p.join()

	best_new_gene = ""
	best_new_num_correct = 0
	for run_index in range(0,len(gene_search_list),num_per_process):
		json_name = "gene_list_len_" + str(len(gene_best_so_far_list)) + "_run_index_" + str(run_index) + ".json"
		with open(json_name, "r") as read_file:
			new_genes_dict = json.load(read_file)
			ListOfSortedGeneTuples = sorted(new_genes_dict.items(), reverse=True, key=lambda v : v[1])	
			if (ListOfSortedGeneTuples[0][1] > best_new_num_correct):
				best_new_num_correct = ListOfSortedGeneTuples[0][1]	
				best_new_gene = ListOfSortedGeneTuples[0][0]

	print("Level " + str(len(gene_best_so_far_list)) + " is done")
	print(best_new_gene)
	print(best_new_num_correct)
	if (best_new_num_correct > gene_best_num_correct):
		gene_best_num_correct = best_new_num_correct
		gene_best_so_far_list.append(best_new_gene)
		gene_search_list.remove(best_new_gene)
		print("Continuing the feature search with this list")
		print(gene_best_so_far_list)
	else:
		print("This iteration did not improve the accuracy, done with search")
		print(gene_best_so_far_list)
		break

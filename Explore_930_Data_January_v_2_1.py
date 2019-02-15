import pandas as pd
import csv
import sys

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

total_tested_dict = dict()
total_correct_dict = dict()

for m in [1, 2, 3, 5, 7, 10, 20, 30]:
	total_tested_dict[m] = 0
	total_correct_dict[m] = 0

for cell_test_name in SC_df_z_test: # go through all the test cells, compare each seperately to the "train" cells
	# get the list of distances for each of the cells
	dist_values_list = []
	print(cell_test_name)
	sys.stdout.flush()
	for cell_train_name in SC_df_z_train:
		arr_1 = SC_df_z_test.loc[:,cell_test_name].values
		arr_2 = SC_df_z_train.loc[:,cell_train_name].values
		dist = euclidean(arr_1, arr_2)
		dist_values_list.append((cell_train_name, dist))
	# now sort them based on distance
	dist_values_list.sort(key=lambda x: x[1])
	# now get accuracies for k-nearest neighbor
	for m in [1, 2, 3, 5, 7, 10, 20, 30]:
		total_G2M = 0
		total_G1 = 0
		total_S = 0
		for k in range(0,m):
			k_cell_train_name = dist_values_list[k][0]
			distance = dist_values_list[k][1]
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
			total_correct_dict[m] += 1
		total_tested_dict[m] += 1

for m in [1, 2, 3, 5, 7, 10, 20, 30]:
	print("Accuracy for " + str(m) + " is " + str(total_correct_dict[m]/total_tested_dict[m]))

import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json


gene_list = list() 
with open('SC_df_z_train_genes.csv', 'r') as f:
    reader = csv.reader(f)
    gene_list = list(reader)[0]
#print(gene_list)

best_features_dict = dict()
with open("cpp_single_features_accuracy_on_training_set.json", "r") as read_file:
	best_features_dict = json.load(read_file)
ListOfSortedGeneTuples = sorted(best_features_dict.items(), reverse=True, key=lambda v : v[1])
print(ListOfSortedGeneTuples)

for k in (1, 3, 5, 7, 10, 20, 30):
	processes = []
	for i in range(0, 252):
		args = ['./Feature_Selection_1-252_Based_On_Single_Order', str(k)]
		for j in range(0, i+1):
			#print(ListOfSortedGeneTuples[j][0]);
			#print(ListOfSortedGeneTuples[j][1]);
			#print(gene_list.index(ListOfSortedGeneTuples[j][0]));
			#print("\n");
			args.append(str(gene_list.index(ListOfSortedGeneTuples[j][0])))
		p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)
		processes.append(p)

	total_correct_list = list()
	for p in processes:
		stdout, stderr = p.communicate()
		print(stderr.decode("utf-8"))
		total_correct_list.append(stdout.decode("utf-8"))

	print(total_correct_list)

	file_name = 'cpp_k_' + str(k) + '_feature_selection_1-252-Based_On_Single_Order_output.csv'
	with open(file_name, 'w', newline='') as myfile:
		wr = csv.writer(myfile) #, quoting=csv.QUOTE_ALL)
		wr.writerow(total_correct_list)
print("finished")

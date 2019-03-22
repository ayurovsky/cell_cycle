import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json
import sys

test_file_prefix_list = ["big_430", "stahlberg", "mESC", "microarray_leg", "bulk_mESC", "big_930"]
#test_file_prefix_list = ["stahlberg"]

train_file_prefix = "microarray_leg_test"

# best features from Recursive feature selection for k=10
features_list = ['TROAP', 'GATA4', 'SERTAD3', 'BCL7A', 'UNG', 'HPS4', 'SLCO1B3', 'DEFB103A', 'HIST1H2AM', 'KRTAP19-1', 'HIST1H4D', 'HIST1H4C', 'HIST1H1A', 'MAGEL2', 'HIST1H4B', 'HIST1H2AJ', 'LCN1', 'HIST1H2AD', 'FGF21', 'FGF20', 'CKS1B', 'HIST1H2AH', 'OR4N4', 'SPP2', 'PTTG2', 'HIST1H4L', 'TP53TG3', 'TRIM43', 'PLCXD1', 'HIST1H3B', 'CTNNBIP1', 'ZNF202', 'EHD1', 'DYRK2', 'BTN3A3', 'EMILIN2', 'MAFK', 'BARD1', 'ARPC5L', 'FBXO22', 'ZNF589', 'HIPK3', 'TNFAIP1', 'TOR3A', 'TAF7', 'CDC40', 'PTTG1IP', 'PAQR4', 'LIMK2', 'CALD1', 'PSAP', 'KIF2C', 'TRIM52', 'SPON2', 'TSPYL4', 'SDC1', 'TDP1', 'VRK1', 'PRPF18', 'TENC1', 'ABHD4', 'RCOR1', 'BFAR', 'KCTD9', 'MTM1', 'PPARGC1B', 'ABTB1', 'EYA3', 'ZDHHC23', 'CD200', 'GNE', 'BIRC3', 'THAP7', 'MTMR4', 'FN3KRP', 'ARHGEF19', 'YWHAZ', 'FBXO32']

canonical_gene_list = list() 
with open(train_file_prefix + '_genes.csv', 'r') as f:
	reader = csv.reader(f)
	canonical_gene_list = list(reader)[0]

canonical_gene_list_indeces = list()
for gene_name in features_list:
   canonical_gene_list_indeces.append(str(canonical_gene_list.index(gene_name)))
#print(canonical_gene_list_indeces)

# when using all features

k = 10 
accuracies = list()
f1_scores = list()
for test_file_prefix in test_file_prefix_list:
	print(test_file_prefix)
	args = ['./Precise_Evaluate_Selected_Features_On_Test', test_file_prefix, train_file_prefix, str(k)]
	args += canonical_gene_list_indeces 
	p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)

	stdout, stderr = p.communicate()
	vals = stdout.decode("utf-8").split("\n")
	err = stderr.decode("utf-8")
	print(err)
	sys.stdout.flush()
	
	accuracies.append(vals[0])
	f1_scores.append(vals[1])	

accs = "\"kNN_MA k=10 RS\",\"" + "\",\"".join(accuracies) + "\""
f1s = "\"kNN_MA k=10 RS\",\"" + "\",\"".join(f1_scores) + "\""
print("Accuracies")
print(accs)
print("F1 Scores")
print(f1s)

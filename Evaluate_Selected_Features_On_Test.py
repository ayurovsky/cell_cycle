import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json


canonical_gene_list = list() 
with open('SC_df_z_train_genes.csv', 'r') as f:
	reader = csv.reader(f)
	canonical_gene_list = list(reader)[0]

features_list = ['HIST2H4A', 'TOP2A', 'HIF1A', 'TLR4', 'CDC25B', 'NR3C1', 'CFB', 'RHNO1', 'IFNA1', 'MCM4', 'IFNB1', 'CEP55', 'C1orf63', 'MYC', 'MAPKAPK2', 'CLTC', 'HJURP', 'TUBB', 'CDC6', 'G2E3', 'MAFG', 'USP1', 'UNG', 'ALAS1', 'VANGL1', 'MAPKAPK5', 'C4A', 'IL1B', 'HMGN1', 'MAPK1']
canonical_gene_list_indeces = list()
for gene_name in features_list:
	canonical_gene_list_indeces.append(str(canonical_gene_list.index(gene_name)))
print(canonical_gene_list_indeces)

k = 10 
args = ['./Evaluate_Selected_Features_On_Test', str(k)]
args += canonical_gene_list_indeces 
p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)

stdout, stderr = p.communicate()
val = stdout.decode("utf-8")

print(val)
print(int(val)/430)

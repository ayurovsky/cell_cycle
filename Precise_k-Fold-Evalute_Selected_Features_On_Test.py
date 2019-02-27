import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json


test_file_prefix = "microarray_leg_test"

canonical_gene_list_indeces = list()
for i in range(0, 1230): #252):
	canonical_gene_list_indeces.append(str(i))	
#print(canonical_gene_list_indeces)

k = 3 
args = ['./Precise_k-Fold-Evaluate_Selected_Features_On_Test', test_file_prefix, str(k)]
args += canonical_gene_list_indeces 
p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)

stdout, stderr = p.communicate()
val = stdout.decode("utf-8")
err = stderr.decode("utf-8")

print(val)
print(err)

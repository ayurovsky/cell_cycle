import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json


phase = "s"
train_file_prefix = "train_precise_" + phase
test_file_prefix = "test_precise"

canonical_gene_list_indeces = list()
for i in range(0, 1230): #252):
	canonical_gene_list_indeces.append(str(i))	
#print(canonical_gene_list_indeces)

k = 5 
total_range = 0.0
if (phase == "s"):
	total_range = 2.0 + (6.283 - 5.9)
elif (phase == "g1"):
	total_range = 5.9 - 4.8 
elif (phase == "g2"):
	total_range = 4.8 - 2
args = ['./Precise_Evaluate_Selected_Features_On_Test', test_file_prefix, train_file_prefix, phase, str(total_range), str(k)]
args += canonical_gene_list_indeces 
p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)

stdout, stderr = p.communicate()
val = stdout.decode("utf-8")
err = stderr.decode("utf-8")

print(val)
print(err)

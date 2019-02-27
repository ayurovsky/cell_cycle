import subprocess
import os
from subprocess import Popen, PIPE
import csv
import json

#test_file_prefix = "bulk_mESC" 
#test_file_prefix = "microarray_leg";
test_file_prefix = "mESC";
#test_file_prefix = "stahlberg";
#test_file_prefix = "big_430";

train_file_prefix = "microarray_leg_test"
#train_file_prefix = "SC_df_z_train"

#canonical_gene_list = list() 
#with open('SC_df_z_train_genes.csv', 'r') as f:
#	reader = csv.reader(f)
#	canonical_gene_list = list(reader)[0]
# all features, sanitiy check 
#features_list = ['HIST2H4A', 'TOP2A', 'HIF1A', 'TLR4', 'CDC25B', 'NR3C1', 'CFB', 'RHNO1', 'IFNA1', 'MCM4', 'IFNB1', 'CEP55', 'C1orf63', 'MYC', 'MAPKAPK2', 'CLTC', 'HJURP', 'TUBB', 'CDC6', 'G2E3', 'MAFG', 'USP1', 'UNG', 'ALAS1', 'VANGL1', 'MAPKAPK5', 'C4A', 'IL1B', 'HMGN1', 'MAPK1', 'G6PD', 'POLR1B', 'KEAP1', 'CDCA8', 'FAM72B', 'CRP', 'LIPH', 'RAC1', 'CCR3', 'CD40LG', 'CKAP5', 'TFRC', 'C3AR1', 'POLR2A', 'ABCF1', 'MAPK3', 'SREK1', 'GNAS', 'RAD21', 'CXCL9', 'HIST2H2BE', 'C15orf23', 'GRPEL1', 'IL6R', 'HLA-DRA', 'LAMC1', 'CXCL5', 'CENPF', 'CCNF', 'CXCR2', 'RIPK1', 'FZR1', 'MKI67', 'CCL16', 'MAP2K1', 'C1R', 'CCR1', 'SLBP', 'PIF1', 'CXCR1', 'CFL1', 'DMXL2', 'CHAF1B', 'CCL11', 'FAM83D', 'IL1RAP', 'IL22', 'GPR37', 'NFKB1', 'HIST1H2AC', 'POLQ', 'MAPK8', 'IL2', 'ORC1', 'CASP8AP2', 'CDC42', 'BUB1', 'ATF2', 'BMP1', 'KIFC1', 'CDKN3', 'IQGAP3', 'RPL19', 'CDCA7', 'MKNK1', 'EEF1G', 'KPNA2', 'ANLN', 'CDCA5', 'CCL24', 'IFNG', 'CKAP2L', 'STAT1', 'MBL2', 'DEPDC1B', 'CKS2', 'OAZ1', 'PPIA', 'RHOA', 'HPRT1', 'EZH2', 'DEPDC1', 'LBR', 'BIRC5', 'CDCA2', 'MAX', 'ARL6IP1', 'C1S', 'NFATC3', 'UBE2C', 'LTB', 'FANCD2', 'GAPDH', 'CDCA3', 'FAM84B', 'ZC3HC1', 'GUSB', 'FAM111B', 'LTA', 'MAPK14', 'INSIG2', 'PCNA', 'CCR7', 'CCL22', 'ROCK2', 'CXCR4', 'PRR11', 'CCNA2', 'RFC4', 'PTTG1', 'NUF2', 'DAXX', 'CASP3', 'IL12B', 'TGFB1', 'CENPQ', 'NUSAP1', 'E2F1', 'VPS25', 'NCAPH', 'TRAF2', 'PSRC1', 'DDX11', 'TNF', 'PGK1', 'C1QB', 'BTBD3', 'FAM189B', 'CCL13', 'ITGB2', 'DLGAP5', 'CCL8', 'RAPGEF2', 'NFE2L2', 'CD55', 'PTK2', 'CDKN1B', 'MCM5', 'IL23A', 'RAF1', 'KIF23', 'PLA2G4A', 'CDC20', 'PLAT', 'MAFF', 'SLC17A2', 'CFD', 'CCL19', 'MEF2BNB-MEF2B', 'ESCO2', 'SMTN', 'IL11', 'H2AFX', 'IL13', 'MEF2D', 'ZNF367', 'CCNB2', 'CXCL1', 'RIPK2', 'CREB1', 'CXCL2', 'SERPING1', 'TPX2', 'JUN', 'MEF2A', 'HSPB2', 'CD40', 'TBP', 'MYD88', 'NEIL3', 'SAPCD2', 'MAP3K7', 'HRAS', 'CCNE1', 'ERBB2', 'LY96', 'VCAM1', 'IL18', 'WDR76', 'CSF2', 'TGFB3', 'HSPB1', 'GMNN', 'MAP2K4', 'MIS18BP1', 'BCL6', 'C5', 'GOLGA8B', 'MBOAT1', 'MSL1', 'MAP3K5', 'PDGFA', 'RPS6KA5', 'DDIT3', 'C6', 'MAP3K1', 'MAP2K6', 'FOS', 'TGFBR1', 'HDAC4', 'IL18RAP', 'MSH5', 'CXCL6', 'TLR1', 'DDX58', 'IL7', 'MAP3K9', 'ITGA1', 'TOLLIP', 'PCF11', 'MEF2C', 'CEBPB', 'TLR3', 'PRKCA', 'GNAQ', 'TTK', 'E2F2', 'TTC38', 'CSF1', 'LIMK1', 'NOX1', 'PLCB1']
# 100 features
#features_list = ['HIST2H4A', 'TOP2A', 'HIF1A', 'TLR4', 'CDC25B', 'NR3C1', 'CFB', 'RHNO1', 'IFNA1', 'MCM4', 'IFNB1', 'CEP55', 'C1orf63', 'MYC', 'MAPKAPK2', 'CLTC', 'HJURP', 'TUBB', 'CDC6', 'G2E3', 'MAFG', 'USP1', 'UNG', 'ALAS1', 'VANGL1', 'MAPKAPK5', 'C4A', 'IL1B', 'HMGN1', 'MAPK1', 'G6PD', 'POLR1B', 'KEAP1', 'CDCA8', 'FAM72B', 'CRP', 'LIPH', 'RAC1', 'CCR3', 'CD40LG', 'CKAP5', 'TFRC', 'C3AR1', 'POLR2A', 'ABCF1', 'MAPK3', 'SREK1', 'GNAS', 'RAD21', 'CXCL9', 'HIST2H2BE', 'C15orf23', 'GRPEL1', 'IL6R', 'HLA-DRA', 'LAMC1', 'CXCL5', 'CENPF', 'CCNF', 'CXCR2', 'RIPK1', 'FZR1', 'MKI67', 'CCL16', 'MAP2K1', 'C1R', 'CCR1', 'SLBP', 'PIF1', 'CXCR1', 'CFL1', 'DMXL2', 'CHAF1B', 'CCL11', 'FAM83D', 'IL1RAP', 'IL22', 'GPR37', 'NFKB1', 'HIST1H2AC', 'POLQ', 'MAPK8', 'IL2', 'ORC1', 'CASP8AP2', 'CDC42', 'BUB1', 'ATF2', 'BMP1', 'KIFC1', 'CDKN3', 'IQGAP3', 'RPL19', 'CDCA7', 'MKNK1', 'EEF1G', 'KPNA2', 'ANLN', 'CDCA5', 'CCL24']
# 50 features_list = ['HIST2H4A', 'TOP2A', 'HIF1A', 'TLR4', 'CDC25B', 'NR3C1', 'CFB', 'RHNO1', 'IFNA1', 'MCM4', 'IFNB1', 'CEP55', 'C1orf63', 'MYC', 'MAPKAPK2', 'CLTC', 'HJURP', 'TUBB', 'CDC6', 'G2E3', 'MAFG', 'USP1', 'UNG', 'ALAS1', 'VANGL1', 'MAPKAPK5', 'C4A', 'IL1B', 'HMGN1', 'MAPK1', 'G6PD', 'POLR1B', 'KEAP1', 'CDCA8', 'FAM72B', 'CRP', 'LIPH', 'RAC1', 'CCR3', 'CD40LG', 'CKAP5', 'TFRC', 'C3AR1', 'POLR2A', 'ABCF1', 'MAPK3', 'SREK1', 'GNAS', 'RAD21', 'CXCL9']
# 30 features_list = ['HIST2H4A', 'TOP2A', 'HIF1A', 'TLR4', 'CDC25B', 'NR3C1', 'CFB', 'RHNO1', 'IFNA1', 'MCM4', 'IFNB1', 'CEP55', 'C1orf63', 'MYC', 'MAPKAPK2', 'CLTC', 'HJURP', 'TUBB', 'CDC6', 'G2E3', 'MAFG', 'USP1', 'UNG', 'ALAS1', 'VANGL1', 'MAPKAPK5', 'C4A', 'IL1B', 'HMGN1', 'MAPK1']
#canonical_gene_list_indeces = list()
#for gene_name in features_list:
#	canonical_gene_list_indeces.append(str(canonical_gene_list.index(gene_name)))
#print(canonical_gene_list_indeces)

canonical_gene_list_indeces = list()
for i in range(0, 1230): #252):
	canonical_gene_list_indeces.append(str(i))	
print(canonical_gene_list_indeces)

k = 20 
args = ['./Generic_Evaluate_Selected_Features_On_Test', test_file_prefix, train_file_prefix, str(k)]
args += canonical_gene_list_indeces 
p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)

stdout, stderr = p.communicate()
val = stdout.decode("utf-8")
err = stderr.decode("utf-8")

print(val)
print(err)

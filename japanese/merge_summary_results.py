import sys
import os
import csv
from collections import OrderedDict 
from pathlib import Path

# 51 genie phenotypes have matching phenotypes in japanese

input_match_manifest_file = os.path.join("..", "manifest_files", "Matching_phe_genie_ukbb_jap.txt")
input_jap_manifest_file = os.path.join("..", "manifest_files", "jap_riken_manifest.txt")
input_genie_result_dirs = ["/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/genie_heritability/convert_plato_output/continous/out", 
						   "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/genie_heritability/convert_plato_output/discrete/out"]
input_jap_gwas_results_dir = "../../genie_download_ukbb_jap_results/jap_results"
out_dir = sys.argv[1]
run_index = int(sys.argv[2]) - 1

SIG_PVALUE = 2
MATCH_ALLELES = True

def flip_allele(a):
	if a == 'A':
		return 'T'
	if a == 'T':
		return 'A'
	if a == 'G':
		return 'C'
	if a == 'C':
		return 'G'

def compare_alleles(a1, a2, b1, b2):
	if a1 == b1 and a2 == b2:
		return True
	if a1 == b2 and a2 == b1:
		return True
	if a1 == flip_allele(b1) and a2 == flip_allele(b2):
		return True
	if a1 == flip_allele(b2) and a2 == flip_allele(b1):
		return True
	return False

# Read japanese manifest file and create a map with Phenotype desc and filename on disk
# The result files were downloaded using the ID number as dir name
jap_phe_to_file_map = {}
with open(input_jap_manifest_file) as inf:
	for line in inf:
		line = line.strip()
		row = line.split('\t')
		row[0] = row[0].strip()
		row[1] = row[1].strip('\"')
		row[1] = row[1].strip()
		in_dir = os.path.join(input_jap_gwas_results_dir, row[0])
		if row[8] == "Yes":
			in_dir = os.path.join(input_jap_gwas_results_dir, "QTL_" + row[0])
		input_files = Path(in_dir).rglob('*.txt')
		for input_file in input_files:
			filename = os.path.basename(input_file)
			if filename != "README.txt" and "chrX" not in filename and "chrx" not in filename:
				if row[1] in jap_phe_to_file_map:
					print("ERROR: Multiple files found for " + row[1] + " in directory " + in_dir + ". Cannot resolve, delete unnecessary files from directory")
					sys.exit(1)
				jap_phe_to_file_map[row[1]] = input_file


# Read the matching jap phenptypes for genie phenotypes manually curated to a map
genie_jap_phe_map = OrderedDict()
with open(input_match_manifest_file) as inf:
	reader = csv.reader(inf, delimiter = '\t')
	next(reader)
	for parts in reader:
		for i in range(1,3):
			parts[i] = parts[i].strip('\"')
			parts[i] = parts[i].strip()
			if parts[i] and "chrX" not in parts[i]:
				if parts[0] not in genie_jap_phe_map:
					genie_jap_phe_map[parts[0]] = set()
				genie_jap_phe_map[parts[0]].add(parts[i])

				if parts[i] not in jap_phe_to_file_map:
					print("Error could not find jap results files for jap phe:" + parts[i])
					sys.exit(1)

# pick a single phe based on index passed as commandline arg 
# this is for distributed runs - run each phe on seperate core
genie_phe = list(genie_jap_phe_map.keys())[run_index]
genie_jap_phe = genie_jap_phe_map[genie_phe]

genie_header_parse_list = ["Outcome", "Var1_ID", "chr", "pos", "a1", "a2", "Var1_MAF", "Num_NonMissing", "Var1_Pval", "Var1_beta", "rsq"]

# parse genie result file and store in hashmap using chrpos as key
genie_res_map = {}

genie_input_file = None
for genie_dir in input_genie_result_dirs:
	if os.path.exists(os.path.join(genie_dir, "result_" + genie_phe)):
		genie_input_file = os.path.join(genie_dir, "result_" + genie_phe)
		break

with open(genie_input_file) as inf:
	reader = csv.reader(inf, delimiter = '\t')
	header_genie = next(reader)
	header_index = []
	for i in range(0, len(genie_header_parse_list)):
		if genie_header_parse_list[i] in header_genie:
			header_index.append(i)
		else:
			header_index.append(None)
	for row in reader:
		if row[0] == genie_phe:
			res_row = []
			for i in header_index:
				if i != None:
					res_row.append(row[i])
				else:
					res_row.append("")
			res_row[2] = res_row[2].replace("chr", "")
			genie_res_map[res_row[2] + ':' + res_row[3]] = res_row


# parse jap result file and store in hashmap using chrpos as key
table_header = [set(["SNP","SNPID","MARKER","Variants"]), set(["CHR"]), set(["POS","BP"]), set(["A1","Allele1","NON_EFFECT_ALLELE","REF"]),
				set(["A2","Allele2","EFFECT_ALLELE","ALT"]), set(["Frq","EFFECT_ALLELE_FREQ","A2Frq","FREQ1","MAF"]), 
				set(["CASE_FREQ1","AF.Cases","Case_frq"]), set(["CTRL_FREQ1","AF.Controls","Ctrl_frq"]), set(["N"]),
				set(["OR"]), set(["BETA"]), set(["P","p.value","PVALUE","P_BOLT"]), set(["Rsq","RSQR"])]


jap_res_map = {}
for jap_phe in genie_jap_phe:
	jap_phe_result_file = jap_phe_to_file_map[jap_phe]
	header_order = []
	with open(jap_phe_result_file) as inf:
		header = next(inf)
		header = header.strip()
		header = header.split()
		for t_header_set in table_header:
			header_not_found = True
			for t_header in t_header_set:
				if t_header in header:
					header_order.append(header.index(t_header))
					header_not_found = False
			if header_not_found:
				header_order.append(None)

		for line in inf:
			line = line.strip()
			parts = line.split()
			row_insert = [jap_phe]
			for i in header_order:
				if i is not None and parts[i] != '.' and parts[i] != '-' and parts[i] != 'NA':
					row_insert.append(parts[i])
				else:
					row_insert.append(None)

			if row_insert[12] and float(row_insert[12]) < SIG_PVALUE:
				row_insert[2] = row_insert[2].replace("chr", "")
				if row_insert[6]:
					row_insert[6] = float(row_insert[6])
					if row_insert[6] > 0.5:
						row_insert[6] = 1 - row_insert[6]
				else:
					row_insert[6] = ""

				if not row_insert[9]:
					row_insert[9] = ""

				chrpos = row_insert[2] + ':' + row_insert[3]
				if chrpos not in jap_res_map or float(jap_res_map[chrpos][5]) > float(row_insert[12]):
					jap_res_map[chrpos] = [row_insert[0], row_insert[4], row_insert[5], 
								str(row_insert[6]), row_insert[9], row_insert[12]]


if not os.path.exists(out_dir):
	os.mkdir(out_dir)

# write the output file
with open(os.path.join(out_dir, genie_phe), 'w') as of:
	of.write('\t'.join(genie_header_parse_list) + '\t' + "jap_phe\tjap_A1\tjap_A2\tjap_MAF\tjap_N\tjap_P\n")
	for chrpos in genie_res_map:
		if chrpos in jap_res_map:
			if MATCH_ALLELES:
				if compare_alleles(genie_res_map[chrpos][4], genie_res_map[chrpos][5], jap_res_map[chrpos][1], jap_res_map[chrpos][2]):
					of.write('\t'.join(genie_res_map[chrpos]) + '\t' + '\t'.join(jap_res_map[chrpos]) + '\n')
				else:
					print('\t'.join(genie_res_map[chrpos]) + '\t' + '\t'.join(jap_res_map[chrpos]) + '\n')
			else:
				of.write('\t'.join(genie_res_map[chrpos]) + '\t' + '\t'.join(jap_res_map[chrpos]) + '\n')




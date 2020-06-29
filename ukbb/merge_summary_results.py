import sys
import os
import csv
from collections import OrderedDict 
from pathlib import Path

# 51 genie phenotypes have matching phenotypes in ukbb

input_match_manifest_file = os.path.join("..", "manifest_files", "Matching_phe_genie_ukbb_jap.txt")
input_ukbb_manifest_file = os.path.join("..", "manifest_files", "UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Manifest 201807.txt")
input_genie_result_dirs = ["/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/genie_heritability/convert_plato_output/continous/out", 
						   "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/genie_heritability/convert_plato_output/discrete/out"]
input_ukbb_gwas_results_dir = "../../genie_download_ukbb_jap_results/ukbb_results"
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

# Read ukbb manifest file and create a map with Phenotype desc and filename on disk
# The result files were downloaded using the ID number as dir name
ukbb_phe_to_file_map = {}
with open(input_ukbb_manifest_file) as inf:
	for i in range(0, 26):
		next(inf)
	for line in inf:
		line = line.strip()
		row = line.split('\t')
		row[0] = row[0].strip()
		row[1] = row[1].strip('\"')
		row[1] = row[1].strip()
		if row[3] == "both_sexes" and not row[0].endswith("irnt"):
			in_dir = os.path.join(input_ukbb_gwas_results_dir, row[0])
			if os.path.exists(in_dir):
				input_files = []
				for in_f in Path(in_dir).rglob('*.tsv'):
					input_files.append(in_f)
				if len(input_files) > 1:
					print("ERROR: Multiple files found for " + row[1] + " in directory " + in_dir + ". Cannot resolve, delete unnecessary files from directory")
					sys.exit(1)
				ukbb_phe_to_file_map[row[1]] = input_files[0]

# Read the matching ukbb phenptypes for genie phenotypes manually curated to a map
genie_ukbb_phe_map = OrderedDict()
with open(input_match_manifest_file) as inf:
	reader = csv.reader(inf, delimiter = '\t')
	next(reader)
	for parts in reader:
		for i in range(3,8):
			parts[i] = parts[i].strip()
			if parts[i]:
				if parts[0] not in genie_ukbb_phe_map:
					genie_ukbb_phe_map[parts[0]] = set()
				genie_ukbb_phe_map[parts[0]].add(parts[i])

				if parts[i] not in ukbb_phe_to_file_map:
					print("Error could not find ukbb results files for ukbb phe:" + parts[i])
					sys.exit(1)

print("total" + str(len(genie_ukbb_phe_map)))


# pick a single phe based on index passed as commandline arg 
# this is for distributed runs - run each phe on seperate core
genie_phe = list(genie_ukbb_phe_map.keys())[run_index]
genie_ukbb_phes = genie_ukbb_phe_map[genie_phe]

# parse genie result file and store in hashmap using chrpos as key
genie_header_parse_list = ["Outcome", "Var1_ID", "chr", "pos", "a1", "a2", "Var1_MAF", "Num_NonMissing", "Var1_Pval", "Var1_beta", "rsq"]

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


# parse ukbb result file and store in hashmap using chrpos as key
ukbb_res_map = {}
for ukbb_phe in genie_ukbb_phes:
	ukbb_phe_result_file = ukbb_phe_to_file_map[ukbb_phe]
	header_order = []
	with open(ukbb_phe_result_file) as inf:
		header = next(inf)
		header = header.strip()
		header = header.split('\t')
		index_remove = []
		if "expected_min_category_minor_AC" in header:
			index_remove.append(header.index("expected_min_category_minor_AC"))
		if "expected_case_minor_AC" in header:
			index_remove.append(header.index("expected_case_minor_AC"))
		index_remove.sort(reverse=True)

		for line in inf:
			line = line.strip()
			parts = line.split('\t')
			for i in index_remove:
				del parts[i]

			# if the variant is not low_confidence_variant
			if parts[3] == "false" and (float(parts[10]) < SIG_PVALUE):
					var_id = parts[0].split(':')
					major_allele = set([var_id[2], var_id[3]])
					major_allele.remove(parts[1])

					chrpos = var_id[0] + ':' + var_id[1]
					if chrpos not in ukbb_res_map or float(ukbb_res_map[chrpos][5]) > float(parts[10]):
						ukbb_res_map[chrpos] = [ukbb_phe, parts[1], major_allele.pop(), parts[2], parts[4], parts[10]]

if not os.path.exists(out_dir):
	os.mkdir(out_dir)

# write the output file
with open(os.path.join(out_dir, genie_phe), 'w') as of:
	of.write('\t'.join(genie_header_parse_list) + '\t' + "ukbb_phe\tukbb_A1\tukbb_A2\tukbb_MAF\tukbb_N\tukbb_P\n")
	for chrpos in genie_res_map:
		if chrpos in ukbb_res_map:
			if MATCH_ALLELES:
				if compare_alleles(genie_res_map[chrpos][4], genie_res_map[chrpos][5], ukbb_res_map[chrpos][1], ukbb_res_map[chrpos][2]):
					of.write('\t'.join(genie_res_map[chrpos]) + '\t' + '\t'.join(ukbb_res_map[chrpos]) + '\n')
				else:
					print('\t'.join(genie_res_map[chrpos]) + '\t' + '\t'.join(ukbb_res_map[chrpos]) + '\n')
			else:
				of.write('\t'.join(genie_res_map[chrpos]) + '\t' + '\t'.join(ukbb_res_map[chrpos]) + '\n')




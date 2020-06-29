import sys
import os
from collections import OrderedDict

genie_jap_input_dir = "../japanese/out_merge_alleles"
genie_ukbb_input_dir = "../ukbb/out_merge_alleles"
genie_jap_ukbb_input_dir = "../japanese_ukbb/out_merge_alleles"
results_chrpos_allele_gene_file = "../pos_gene_annotation/all_positions_alleles_genes"
out_dir = sys.argv[1]
run_index = int(sys.argv[2]) - 1

PVAL_CUTOFF = 0.0001

if not os.path.exists(out_dir):
	os.mkdir(out_dir)

phes_list = ["agitation", "alcohol_heavy", "ALP", "aorta_dilatation", "Appetite_change_increase", "aPTT", "BASOPHIL", "BMI", "brain_aneurysmm", 
		"brain_atherosclerosis", "brain_stenosis", "breast_cancer", "BUN", "cal_score_cont", "CALCIUM", "cataract", "cholecystitis", "CL", "cognitive_funx", 
		"coronary_plaque", "coronary_stenosis", "CREATININE", "daily_coffee01_23", "Depressed_mood", "depression_score", "DEXA_bone_density", "DM_diagnosis", 
		"duodenal_ulcer", "education", "EKG_rate", "EOSINOPHIL", "exercise_continuous", "fatigue", "fatty_liver_4grade0_123", "FEV1_L", "FEV1_percent", "free_T4", 
		"FVC_L", "gastric_cancer", "gastric_ulcer", "GB_stone", "GERD", "GFR", "GGT", "glucose", "GOT", "GPT", "guilt_feeling", "HB", "HBA1C", "HBV", "HCT", 
		"HCV", "HDL", "Height", "HTN_diagnosis", "Inbody_Body_fat_mass", "Inbody_Fat_percent", "IOP_left", "IOP_rt", "K", "LDL_chol", "loss_of_interest", "LYMPHOCYTE", 
		"macular_change", "MCH", "MCHC", "MCV", "metabolic_syndrome", "MONOCYTE", "MPV", "Na", "optic_fiber_loss", "PCT", "PHOSPHORUS", "PLT", "PT", "RBC", "RDW", 
		"Renal_stone", "retardation", "SEG_NEUT.", "serum_albumin", "serum_protein", "serum_total_bil", "smoking_2", "SOL", "spine_compression_fracture", 
		"spine_disc_narrowing", "spine_spondylolisthesis", "spine_spondylosis", "suicidal", "TG", "thyroid_cancer", "Total_chol", "total_fat_CT_mm2", "TSH", 
		"URIC_ACID", "urine_PU_cat", "VitD3", "WASO", "WBC", "WC", "Weight", "urine_PH"]


class PvalCols:
	def __init__(self, sig_cols, not_sig_cols):
		self.sig_cols = sig_cols
		self.not_sig_cols = not_sig_cols

	def get_sig_cols(self):
		return self.sig_cols

	def get_not_sig_cols(self):
		return self.not_sig_cols

	def set_sig_cols_indices(self, indices):
		self.sig_cols_indices = indices

	def set_not_sig_cols_indices(self, indices):
		self.not_sig_cols_indices = indices

	def get_sig_cols_indices(self):
		return self.sig_cols_indices

	def get_not_sig_cols_indices(self):
		return self.not_sig_cols_indices


def find_header_indices(header, string_arr):
	indices = []
	for h in string_arr:
		if h in header:
			indices.append(header.index(h))
		else:
			print("Error Could not find index for:" + h + " in " + str(header))
			sys.exit(1)
	return indices

def count_pval_signifcant(row, indices):
	count = 0
	for pval_header_index in indices:
		if float(row[pval_header_index]) < PVAL_CUTOFF:
			count += 1
	return count

def calculate_stats_from_file(in_file, pvalcols_arr, chrpos_gene_map):
	if not os.path.exists(in_file):
		return (None, None)
	with open(in_file) as inf:
		header = next(inf)
		header = header.strip()
		header = header.split('\t')
		chr_index = header.index("chr")
		pos_index = header.index("pos")

		col_index_set = set()
		not_sig_col_set = set()
		for pvalcols_ob in pvalcols_arr:
			pvalcols_ob.set_sig_cols_indices(find_header_indices(header, pvalcols_ob.get_sig_cols()))
			pvalcols_ob.set_not_sig_cols_indices(find_header_indices(header, pvalcols_ob.get_not_sig_cols()))
			col_index_set.update(pvalcols_ob.get_sig_cols_indices())
			col_index_set.update(pvalcols_ob.get_not_sig_cols_indices())

		sig_p_count_list = [0] * (len(pvalcols_arr) + 2)
		col_index_genes_map = {}
		all_genes = set()
		sig_genes_set_rm_list = [set() for i in range(len(pvalcols_arr) + 2)]
		for line in inf:
			line = line.strip()
			row = line.split('\t')
			sig_p_count_list[0] += 1
			genes = chrpos_gene_map[row[chr_index] + ':' + row[pos_index]]
			if genes:
				all_genes.update(genes)
			for i in range(0, len(pvalcols_arr)):
				s_count = count_pval_signifcant(row, pvalcols_arr[i].get_sig_cols_indices())
				ns_count = count_pval_signifcant(row, pvalcols_arr[i].get_not_sig_cols_indices())
				# Count any significant snp and gene
				if i == 0 and (s_count > 0 or ns_count > 0):
					sig_p_count_list[1] += 1
				# Count exclusively significant in the group
				if s_count == len(pvalcols_arr[i].get_sig_cols_indices()) and ns_count == 0:
					sig_p_count_list[i + 2] += 1
			
			for col_index in col_index_set:
				if genes and float(row[col_index]) < PVAL_CUTOFF:
					if col_index not in col_index_genes_map:
						col_index_genes_map[col_index] = set()
					col_index_genes_map[col_index].update(genes)

		sig_genes_set_list = [set() for i in range(len(pvalcols_arr) + 2)]
		sig_genes_set_list[0]= all_genes
		for i in range(0, len(pvalcols_arr)):
			for col_index in set(pvalcols_arr[i].get_sig_cols_indices() + pvalcols_arr[i].get_not_sig_cols_indices()):
				sig_genes_set_list[1].update(col_index_genes_map[col_index])

			intersect_set = None
			for index in pvalcols_arr[i].get_sig_cols_indices():
				if not intersect_set:
					intersect_set = col_index_genes_map[index]
				else:
					intersect_set = intersect_set & col_index_genes_map[index]

			non_sig_gene_set = set()
			for index in pvalcols_arr[i].get_not_sig_cols_indices():
				non_sig_gene_set.update(col_index_genes_map[index])

			sig_genes_set_list[i + 2] = intersect_set - non_sig_gene_set

		return (sig_p_count_list, sig_genes_set_list)


def parse_chrpos_genes(anno_file):
	chrpos_gene = {}
	with open(anno_file) as inf:
		next(inf)
		for line in inf:
			row = line.split('\t')
			row[len(row) - 1] = row[len(row) - 1].strip()
			chrpos = row[0] + ':' + row[1]
			chrpos_gene[chrpos] = set()
			if row[4]:
				for gene in row[4].split(','):
					chrpos_gene[chrpos].add(gene)
	return chrpos_gene


def write_to_file(num_fields, count_list, gene_set_list, gene_fp, stat_fp):
	if count_list == None or gene_set_list == None:
		for i in range(0, num_fields + 2):
			stat_fp.write(",,")
			if i > 1:
				gene_fp.write(',')
	else:
		for i in range(0, len(count_list)):
			stat_fp.write(',' + str(count_list[i]))
			stat_fp.write(',' + str(len(gene_set_list[i])))

			if i > 1:
				gene_fp.write(',"' + '|'.join(gene_set_list[i]) + '"')


chrpos_gene_map = parse_chrpos_genes(results_chrpos_allele_gene_file)

out_genes_file = os.path.join(out_dir, "genes_" + phes_list[run_index] + ".csv")
out_stat_file = os.path.join(out_dir, "stat_" + phes_list[run_index] + ".csv")
phe = phes_list[run_index]

with open(out_genes_file, 'w') as ogf, open(out_stat_file, 'w') as osf:
	
	print("Processing phe:" + phe)
	ogf.write(phe)
	osf.write(phe)
	
	genie_jap_file = os.path.join(genie_jap_input_dir, phe)
	pcols = [PvalCols(["Var1_Pval"], ["jap_P"]), PvalCols(["jap_P"], ["Var1_Pval"]), PvalCols(["Var1_Pval", "jap_P"], [])]
	count_list, gene_set_list = calculate_stats_from_file(genie_jap_file, pcols, chrpos_gene_map)
	write_to_file(len(pcols), count_list, gene_set_list, ogf, osf)

	genie_ukbb_file = os.path.join(genie_ukbb_input_dir, phe)
	pcols = [PvalCols(["Var1_Pval"], ["ukbb_P"]), PvalCols(["ukbb_P"], ["Var1_Pval"]), PvalCols(["Var1_Pval", "ukbb_P"], [])]
	count_list, gene_set_list = calculate_stats_from_file(genie_ukbb_file, pcols, chrpos_gene_map)
	write_to_file(len(pcols), count_list, gene_set_list, ogf, osf)

	genie_jap_ukbb_file = os.path.join(genie_jap_ukbb_input_dir, phe)
	pcols = [PvalCols(["Var1_Pval"], ["ukbb_P", "jap_P"]), PvalCols(["Var1_Pval", "ukbb_P", "jap_P"], [])]
	count_list, gene_set_list = calculate_stats_from_file(genie_jap_ukbb_file, pcols, chrpos_gene_map)
	write_to_file(len(pcols), count_list, gene_set_list, ogf, osf)

	ogf.write('\n')
	osf.write('\n')





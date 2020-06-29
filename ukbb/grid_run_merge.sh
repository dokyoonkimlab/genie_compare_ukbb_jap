#!/bin/bash
#$ -N ukbb_merge
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=1:00:00
#$ -t 1-102

module load Anaconda/3.7
python merge_summary_results.py out_merge_alleles $SGE_TASK_ID



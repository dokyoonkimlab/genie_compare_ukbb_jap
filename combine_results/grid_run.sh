#!/bin/bash
#$ -N stat_combine
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=0:30:00
#$ -t 1-105

python generate_stats.py out $SGE_TASK_ID


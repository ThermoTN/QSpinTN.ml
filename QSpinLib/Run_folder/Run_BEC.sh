#!/bin/bash 
#SBATCH --job-name=Delta=0.5_GSDY-YC6x15-D2000 
#SBATCH --mem=40gb 
#SBATCH --partition=cpu_2020 
#SBATCH --time=99-00:00:00 
#SBATCH --nodes=1 
#SBATCH --exclude=node1078 
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1 
#SBATCH --output=/public/home/yuangao/QSpace_v4.1/tanTRG_KCSO_Th/Run_folder/log_BEC/YC6x12_D=1000_Jxy=0.1_H=2.8.log 
cd /public/home/yuangao/QSpace_v4.1/tanTRG_KCSO_Th/ 
matlab -nodesktop -nodisplay -r "Run" 

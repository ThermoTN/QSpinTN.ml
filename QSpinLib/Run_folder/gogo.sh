#!/bin/bash 
#SBATCH --job-name=gg
#SBATCH --mem=4gb 
#SBATCH --partition=cpu_2020 
#SBATCH --time=14-00:00:00 
#SBATCH --nodes=1 
#SBATCH --exclude=node2029
#SBATCH --ntasks-per-node=4 
#SBATCH --cpus-per-task=1 
#SBATCH --output=/public2/home1/yuangao/QSpinLibQMagen/QSpinLib_DMRG_SCAN/Run_folder/gg.log 
cd /public2/home1/yuangao/QSpinLibQMagen/QSpinLib_DMRG_SCAN/Run_folder
matlab -nodesktop -nodisplay -r "go_script"

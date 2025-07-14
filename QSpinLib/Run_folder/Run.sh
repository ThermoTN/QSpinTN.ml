#!/bin/bash 
#SBATCH --job-name=NoSymm_TriYC6x12 
#SBATCH --mem=80gb 
#SBATCH --partition=cpu_2020 
#SBATCH --time=14-00:00:00 
#SBATCH --nodes=1 
#SBATCH --exclude=node1065 
#SBATCH --ntasks-per-node=4 
#SBATCH --cpus-per-task=1 
#SBATCH --output=/public2/home1/yuangao/QSpinLibQMagen/QSpinLib/Run_folder/log/YC6x12_Tri_D=800.log 
cd /public2/home1/yuangao/QSpinLibQMagen/QSpinLib
matlab -nodesktop -nodisplay -r "RunQSpinLib_TriHei" 

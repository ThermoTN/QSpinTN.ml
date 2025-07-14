% Q = [0.065 0.09 0.115 0.14 0.165 0.19 0.2025 0.215 0.2275 0.24 0.2525 0.265 0.29 0.315 0.34 0.365 0.415]; %[0:2:14 20 32 36 54:2:74];
Q = 0:2:20;
% Q = 25;
% Ac = {'zz', 'xx', 'yy'};
% for ItK = 1:3
for It = 1:length(Q)
    fileID = fopen('Run.sh','w');
    fprintf(fileID,'#!/bin/bash \n');
    fprintf(fileID,'#SBATCH --job-name=Delta=0.5_GSDY-YC6x15-D2000 \n');
    fprintf(fileID,'#SBATCH --mem=40gb \n');
    fprintf(fileID,'#SBATCH --partition=cpu_2020 \n');
    fprintf(fileID,'#SBATCH --time=99-00:00:00 \n');
    fprintf(fileID,'#SBATCH --nodes=1 \n');
    fprintf(fileID,'#SBATCH --exclude=node1078 \n');
    fprintf(fileID,'#SBATCH --ntasks-per-node=8 \n');
    fprintf(fileID,'#SBATCH --cpus-per-task=1 \n');
    % fprintf(fileID,'#SBATCH --output=/public/home/yuangao/QSpace_v4.1/tanTRG_KCSO_Opt/Run_folder/log/Jxy=0.11085-g=8.3948-ES=30.0489-B=%f.log \n', Q(It));
    fprintf(fileID,'#SBATCH --output=/public/home/yuangao/QSpace_v4.1/tanTRG_KCSO_Opt/Run_folder/log/Jxy=0.0705-g=7.8-ES=34.5828-B=%f.log \n', Q(It));
    fprintf(fileID,'cd /public/home/yuangao/QSpace_v4.1/tanTRG_KCSO_Opt/ \n');
    % fprintf(fileID,['matlab -nodesktop -nodisplay -r "CalRslt(0.11085, 8.3948, 30.0489, %f)" \n'], Q(It))
    fprintf(fileID,['matlab -nodesktop -nodisplay -r "CalRslt(0.0705, 7.8, 34.5828, %f)" \n'], Q(It))
    system('sbatch Run.sh')
    pause(2);
end
% end

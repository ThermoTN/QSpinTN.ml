clear; clc;
maxNumCompThreads(4);
addpath(genpath('../../ManyBodySolver/'));
addpath(genpath('../../SpinModel/'));
addpath(genpath('../../svd_lapack_interface/'));
% set model
Para.IntrcMap_Name = 'IntrcMap_XXZtest';

% set coupling strength
Para.Model.Jxy = 1;
Para.Model.Jz = 1;

Para.d = 2; % d = 2S + 1 with S = 1/2

% set lattice geometry
Para.Geo.L = 12;
Para.Geo.BC = 'OBC';
Para.L = Para.Geo.L;

% set magnetic field
Para.Field.h = [0,0,0.1];

%
Para.ThDQ = 'Cm'; % 'Cm': sepcific heat, 'Chi': magnetic susceptibility (M/h) [condition: norm(Para.Field.h) ~= 0]
Para.ManyBodySolver = 'ED'; % 'ED': exact diagonalization, 'ED_SM': construct sparse matrix of Hamiltonian
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

switch Para.ManyBodySolver
    case 'ED'
        Rslt = GetEDRslt(Para);
    case 'ED_SM'
        H = GetEDSMRslt(Para);
end


% Dispaly result
if strcmp(Para.ManyBodySolver, 'ED')
    switch Para.ThDQ
        case 'Cm'
            semilogx(Rslt.T, Rslt.Cm, '-o', 'linewidth', 2);
            xlabel('$T$', 'Interpreter', 'latex')
            ylabel('$C_{\rm m}$', 'Interpreter', 'latex')
            set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)
        case 'Chi'
            semilogx(Rslt.T, Rslt.M/norm(Para.Field.h), '-o', 'linewidth', 2);
            xlabel('$T$', 'Interpreter', 'latex')
            ylabel('$\chi$', 'Interpreter', 'latex')
            set(gca, 'XColor', 'k', 'YColor', 'k', 'fontsize', 20, 'fontname', 'times new roman', 'linewidth', 1.5)
    end
end

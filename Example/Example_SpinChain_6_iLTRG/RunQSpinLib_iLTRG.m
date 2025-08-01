% Note: not recommend
clear; clc;
maxNumCompThreads(4);
addpath(genpath('../../ManyBodySolver/'));
addpath(genpath('../../SpinModel/'));
addpath(genpath('../../svd_lapack_interface/'));
% set model
Para.Trotter_Name = 'Trotter_XXZtest';

% set coupling strength
Para.Model.Jxy = 1;
Para.Model.Jz = 0;

Para.d = 2; % d = 2S + 1 with S = 1/2

% set lattice geometry

% set magnetic field
Para.Field.h = [0,0,0.1];

%
Para.ThDQ = 'Cm&Chi';
Para.WorkingMode = 'ThDQ';
Para.ManyBodySolver = 'iLTRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

Rslt = GetiLTRGRslt(Para, Para.ThDQ);

figure(1)
plot(Rslt.T, Rslt.Cm, '-o', 'linewidth', 2);
xlim([0, 4])
xlabel('$T$', 'Interpreter', 'latex')
ylabel('$C$', 'Interpreter', 'latex')
set(gca, 'fontsize', 22, 'fontname', 'times new roman', 'linewidth', 1.5, ...
    'XColor', 'k', 'YColor', 'k')

if norm(Para.Field.h) ~= 0
    figure(2)
    plot(Rslt.T, Rslt.M/norm(Para.Field.h), '-o', 'linewidth', 2);
    xlim([0, 4])
    xlabel('$T$', 'Interpreter', 'latex')
    ylabel('$\chi$', 'Interpreter', 'latex')
    set(gca, 'fontsize', 22, 'fontname', 'times new roman', 'linewidth', 1.5, ...
        'XColor', 'k', 'YColor', 'k')
end
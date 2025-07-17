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
Para.Model.Jz = 1;

Para.d = 2; % d = 2S + 1 with S = 1/2

% set lattice geometry

% set magnetic field
Para.Field.h = [0,0,0.1];

%
Para.ThDQ = 'Cm';
Para.WorkingMode = 'ThDQ';
Para.ManyBodySolver = 'iLTRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

Rslt = GetiLTRGRslt(Para, Para.ThDQ);


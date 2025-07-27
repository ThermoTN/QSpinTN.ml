clear; clc;
maxNumCompThreads(4);
addpath(genpath('../../ManyBodySolver/'));
addpath(genpath('../../SpinModel/'));
addpath(genpath('../../svd_lapack_interface/'));
% set model
Para.IntrcMap_Name = 'IntrcMap_XXZtest';

% set coupling strength
Para.Model.Jxy = 1;
Para.Model.Jz = 0;

Para.d = 2; % d = 2S + 1 with S = 1/2

% set lattice geometry
Para.Geo.L = 50;
Para.Geo.BC = 'OBC';
Para.L = Para.Geo.L;

% set magnetic field
Para.Field.h = [0,0,0];

Para.saveInfo = ['XTRG'];
%
Para.ThDQ = '';
Para.ManyBodySolver = 'XTRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

Rslt = GetXTRGRslt(Para);
savetanTRGRslt(Para, Rslt)

clear; clc;
maxNumCompThreads(4);
addpath(genpath('../../ManyBodySolver/'));
addpath(genpath('../../SpinModel/'));
addpath(genpath('../../svd_lapack_interface/'));
% set model
Para.IntrcMap_Name = 'IntrcMap_TLARX';
% Para.saveInfo = ['debug'];
% set coupling strength
Para.Model.J1xy = 1;
Para.Model.J1z = 1;
Para.Model.JPD = 0.5;
Para.Model.JGamma = 0;
Para.Model.J2xy = 0;
Para.Model.J2z = 0;


Para.d = 2;

% set lattice geometry
Para.Geo.L = 12;
Para.Geo.Lx = 4;
Para.Geo.Ly = 3;
Para.Geo.BCX = 'OBC';
Para.Geo.BCY = 'PBC';
Para.L = Para.Geo.L;

% set magnetic field
Para.Field.h = [0,0,0.1];

%
Para.ThDQ = 'Cm';
Para.ManyBodySolver = 'tanTRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);

Rslt = GettanTRGRslt(Para);


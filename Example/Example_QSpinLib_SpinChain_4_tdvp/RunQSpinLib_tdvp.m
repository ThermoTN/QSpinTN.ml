clear all
clc
addpath(genpath('../../ManyBodySolver/'));
addpath(genpath('../../SpinModel/'));
addpath(genpath('../../svd_lapack_interface/'));

load MPS/'Model=XXZtest-Jxy=1-Jz=1-h=[0 0 0.1]-L=10-BC=OBC-D=100.mat'
maxNumCompThreads(8);
addpath(genpath('../ManyBodySolver/'));

Para.DY.k = pi;
Para.DY.QulocOp = [-0.5, 0; 0, 0.5];

Para.saveInfo = ['k=', mat2str(Para.DY.k, 3), '-Op=Sz'];

Para.tau = 0.1;
Para.MCrit = 100;
Para.E0 = Rslt.E0;
Para.tmax = 10;
Para.t_switch = 2;
Para.Full2 = false;
Para.Stol = 1e-12;

Para = InitKSpaceIntr(Para);
GettdvpRslt(Para, Rslt.T, Rslt.H )
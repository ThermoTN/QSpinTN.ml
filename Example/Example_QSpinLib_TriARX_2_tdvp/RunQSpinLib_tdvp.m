clear all
clc

load MPS/'Model=TLARX-J1xy=1-J1z=1-JPD=0.5-JGamma=0-J2xy=0-J2z=0-h=[0 0 0.1]-L=12-Lx=4-Ly=3-BCX=OBC-BCY=PBC-D=100.mat'
maxNumCompThreads(4);
addpath(genpath('../ManyBodySolver/'));

Para.DY.k = [4/3, 0] * pi;
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
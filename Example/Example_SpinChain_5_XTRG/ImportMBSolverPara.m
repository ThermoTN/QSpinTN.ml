function [ Para ] = ImportMBSolverPara( Para )

% // SETTEN initial tau
Para.tau0 = 0.00025;

% It: XTRG iterations
Para.It = 20;
% //MCrit: bond dimension compressing H^n
% used in SETTN initialization of rho(tau)
% Recommended value 200
% Bigger is more accurate
Para.MCrit = 100;

% //D: bond dimension of rho(beta)
% D_list = [Di for i steps, Dj for j steps, etc .]
% Recommended value 200
% Bigger is more accurate
Para.D_list = 1:1:Para.It;
Para.D_list(:) = 100;
Para.D_list(end:-1:end-4) = 100;

Para.sweep_step = 50;

% //XTRG runtime parameters
Para.Ver = 'Memory';
% max iterations of MPO varitional product
Para.VariProd_step_max = 10000;
% max iterations of MPO varitional sum
Para.VariSum_step_max = 10000;
% max expensian order of SETTN
Para.SETTN_init_step_max = 10000;

end

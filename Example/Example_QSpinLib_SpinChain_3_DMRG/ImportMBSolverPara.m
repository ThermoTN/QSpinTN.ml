function [ Para ] = ImportMBSolverPara( Para )

Para.MCrit = 100;

% //DMRG runtime parameters
Para.Ver = 'Memory';
% max iterations of MPO varitional product
Para.VariProd_step_max = 10000;
% max iterations of MPO varitional sum
Para.VariSum_step_max = 10000;
% max expensian order of SETTN
Para.SETTN_init_step_max = 10000;

Para.DK = 20;
Para.tol = 1e-8;
Para.DMRGStepMax = 100;

end


function [ Para ] = ImportMBSolverPara( Para )

Para.tau0 = 2^-15;
Para.beta_max = 32;
[Para.tau_step, Para.beta_c ] = GetStepFunc(Para.tau0, [1, 0.25], Para.beta_max);

Para.TSRGStol = 1e-12;
Para.Full2 = true; % Full 2-site tdvp
Para.beta_switch = 2; 

Para.MCrit = 200;

Para.Ver = 'Memory';
% max iterations of MPO varitional product
Para.VariProd_step_max = 10000;
% max iterations of MPO varitional sum
Para.VariSum_step_max = 10000;
% max expensian order of SETTN
Para.SETTN_init_step_max = 4;

end

function [ tau_step, beta_c ] = GetStepFunc(tau, cric, beta_list)
Step = floor(log(cric(1)/tau)/log(2));
beta_c = tau * 2.^(2:Step);
tau_step = tau * 2.^(1:Step-1)/4;
beta_c = beta_c * 0.999;
beta_c = [beta_c, cric(1), beta_list(1:end)];
tau_step = [tau_step, (cric(1)-tau*2.^Step)/4, cric(2:end)];
end

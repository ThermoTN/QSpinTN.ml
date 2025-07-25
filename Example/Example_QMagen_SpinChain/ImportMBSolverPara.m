function [ Para ] = ImportMBSolverPara( Para )
% function [ Para ] = ImportMBSolverPara( Para )
switch Para.ManyBodySolver
    case {'ED', 'ED_C'}
        % // beta list for ED
        beta_list = 0.0025.*2.^(0:1:20).*2^0.1;
        for int = 0.3:0.2:0.9
            beta_list = [beta_list, 0.0025.*2.^(0:1:20).*2^int];
        end
        beta_list = sort(beta_list);
        Para.beta_list = beta_list;
        
    case 'iLTRG'
        % // Trotter step
        Para.tau0 = 0.025;
        
        % Percentage change of field when calculate susceptibility
        % Recommended value
        %       Zero field     0.25
        %       Finite field   0.01
        Para.DeltaHRatio = 0.25;
        % // N_max: iLTRG iterations
        Para.N_max = floor(Para.beta_max / 2 / Para.tau) + 3;
        
        % // Trotter order
        Para.TroOrd = '1';      % '1' only
        % // D_max: bond dimension of rho(beta/2)
        Para.D_max = 40;
        
        % // Number of interp
        Para.InterNum = 80;
        
    case 'XTRG'
        Para.RunParallel = true;
        Para.RunParallelNum = 2;
        Para.RunParallelmaxNumCompThreads = 4;
        % // SETTEN initial tau
        Para.tau0 = 0.00025;
        
        % It: XTRG iterations
        Para.It = floor(log(Para.beta_max/Para.tau0)/log(2));
        
        % //D: bond dimension of rho(beta)
        % D_list = [Di for i steps, Dj for j steps, etc .]
        % Recommended value 200
        % Bigger is more accurate
        Para.D_list = 1:1:Para.It;
        Para.D_list(:) = 50;
        Para.D_list(end:-1:end-4) = 50;
        
        Para.sweep_step = 3;
        % //MCrit: bond dimension compressing H^n
        % used in SETTN initialization of rho(tau)
        % Recommended value 200
        % Bigger is more accurate
        Para.MCrit = 50;
        
        % //XTRG runtime parameters
        Para.Ver = 'Memory';
        % max iterations of MPO varitional product
        Para.VariProd_step_max = 10000;
        % max iterations of MPO varitional sum
        Para.VariSum_step_max = 10000;
        % max expensian order of SETTN
        Para.SETTN_init_step_max = 10000;
        
    case 'TSRG'
        
        % // SETTEN initial tau
        Para.tau0 = 2^-15;
        Para.TSRGStol = 1e-8;
        if ~isfield(Para, 'beta_max')
            Para.beta_max = 1;
        else
            Para.beta_max = Para.beta_max + 8 * 4 * 0.125;
        end
        % Para.max_beta = 1000;
        [Para.tau_step, Para.beta_c ] = GetStepFunc(Para.tau0, [0.5, 0.125], [Para.beta_max]);
        Para.Full2 = false;
        % //D: bond dimension of rho(beta)
        % D_list = [Di for i steps, Dj for j steps, etc .]
        % Recommended value 200
        % Bigger is more accurate
        
        % //MCrit: bond dimension compressing H^n
        % used in SETTN initialization of rho(tau)
        % Recommended value 200
        % Bigger is more accurate
        Para.MCrit = 300;
        
        % //XTRG runtime parameters
        Para.Ver = 'Memory';
        % max iterations of MPO varitional product
        Para.VariProd_step_max = 10000;
        % max iterations of MPO varitional sum
        Para.VariSum_step_max = 10000;
        % max expensian order of SETTN
        Para.SETTN_init_step_max = 2;
end

end

function [ tau_step, beta_c ] = GetStepFunc(tau, cric, beta_list)
Step = floor(log(cric(1)/tau)/log(2));
beta_c = tau * 2.^(2:Step);
tau_step = tau * 2.^(1:Step-1)/4;
beta_c = beta_c * 0.999;
beta_c = [beta_c, cric(1), beta_list(1:end)];
tau_step = [tau_step, (cric(1)-tau*2.^Step)/4, cric(2:end)];
end

function [ Para ] = ImportMBSolverPara( Para )

% // Trotter step
Para.tau = 0.025;
Para.beta_max = 30;
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
Para.D_max = 200;

% // Number of interp
Para.InterNum = 80;

Para.TStr_log = datestr(now,'YYYYmmDD_HHMMSS');
end


function [ Para ] = ImportMBSolverPara( Para )

Para.MCrit = 100;
Para.saveMPS = true;
Para.DK = 20;
Para.tol = 1e-8;
Para.DMRGStepMax = 100; % Maximum sweeps per stage
Para.DMRGSweepTol = 1e-8; % Relative energy AND change of maximum bond entropy
Para.DMRGUpdateSites = "auto"; % "auto": two-site -> one-site on energy rise; or 1 / 2
Para.DMRGEnvStorage = 'memory'; % 'memory' or 'disk', for both stages
Para.DMRGEnvRoot = fullfile(tempdir, 'qspintn_dmrg'); % Disk cache parent
Para.DMRGKeepEnv = true; % Retain disk environments; false cleans them on return

end


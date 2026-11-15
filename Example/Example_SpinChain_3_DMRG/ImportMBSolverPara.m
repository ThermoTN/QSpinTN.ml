function [ Para ] = ImportMBSolverPara( Para )

Para.MCrit = 100;
Para.saveMPS = false; % true includes T and H in Rslt; does not write a file
Para.DK = 20;
Para.tol = 1e-12;
Para.DMRGStepMax = 100; % Maximum sweeps per stage
Para.DMRGSweepTol = 1e-8; % Relative energy AND change of maximum bond entropy
Para.DMRGUpdateSites = "auto"; % "auto": two-site -> one-site on energy rise; or 1 / 2 / "cbe" (RSVD)
Para.CBERatio = 0.1; % RSVD expansion target: (1+ratio)*MCrit
Para.CBEOversampling = 8;
Para.CBEPowerIterations = 1;
Para.CBESvdTol = 1e-12; % Expansion singular-value cutoff; DMRG SVD cutoff
Para.CBEExpansionCap = []; % Optional temporary bond cap, >= MCrit
Para.DMRGEnvStorage = 'memory'; % 'memory' or 'disk', for both stages
Para.DMRGEnvRoot = fullfile(tempdir, 'qspintn_dmrg'); % Disk cache parent
Para.DMRGKeepEnv = true; % Retain disk environments; false cleans them on return
Para.DMRGInitialMPSFile = ''; % Previous MAT filename (Rslt.T or T); empty starts fresh

end


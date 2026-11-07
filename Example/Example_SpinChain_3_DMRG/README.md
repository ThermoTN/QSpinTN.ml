# Example: Spin Chain DMRG

Run `RunQSpinLib_DMRG.m` from this directory for the spin-1/2 XXZ chain:

$$H=\sum_{i=1}^{L-1} J_{xy}(S_i^x S_{i+1}^x+S_i^y S_{i+1}^y)+J_z S_i^z S_{i+1}^z-\sum_{i=1}^L\mathbf{h}\cdot\mathbf{S}_i.$$

## Parameters

`ImportMBSolverPara.m` sets the complete DMRG configuration:

```matlab
Para.MCrit = 100;
Para.saveMPS = false; % true returns T/H; it does not write a file
Para.DK = 20;
Para.tol = 1e-12; % Local Lanczos Ritz-energy change / breakdown
Para.DMRGStepMax = 100; % Maximum sweeps per stage
Para.DMRGSweepTol = 1e-8; % Relative energy AND change of maximum bond entropy
Para.DMRGUpdateSites = "auto"; % Or numeric 1 / 2
Para.DMRGEnvStorage = 'memory'; % Or 'disk'
Para.DMRGEnvRoot = fullfile(tempdir, 'qspintn_dmrg');
Para.DMRGKeepEnv = true; % false removes this run's disk cache on exit
```

Auto starts with two-site sweeps and switches to one-site on a resolved energy
rise. Each stage has its own `DMRGStepMax` budget; convergence can occur before a
switch. To exercise environment saving, select `'disk'`; every call creates a
unique directory, reported in `Rslt.EnvDir`. Environment files alone cannot
resume a run. See the [solver contract](../../ManyBodySolver/DMRG/README.md) for
exact convergence rules, initialization limits, timing fields and cleanup.

## Results

```matlab
Para.ManyBodySolver = 'DMRG';
Para.fileID = 1;
Para = ImportMBSolverPara(Para);
Rslt = GetDMRGRslt(Para);
disp(Rslt.DMRGInfo)
```

- `E0`: total ground-state energy estimate in natural units.
- `M`: per-site magnetization along the field; assigned zero at zero field.
- `EE`: final maximum bond entropy (pre-truncation in two-site mode).
- `DMRGInfo`: convergence/stop reason, sweep counts, switch diagnostics and timings.
- `T`, `H`: included when `saveMPS=true`, for an explicit save and TDVP preparation.

Check `DMRGInfo.Converged` and `StopReason` and assess bond/Krylov convergence.
The script leaves results in memory and does not export figures or files.

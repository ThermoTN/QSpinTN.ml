# Example\: Ground-state dynamical spin structure factor of 1D Heisenberg Chain
We offer here an example of dynamical spin structure factor (DSS). \
**Note: Please read** ```Example_SpinChain_4_tdvp``` **firstly**.\
Run **RunQSpinLib_DMRG.m** firstly to obtain the groud-state MPS. \
Run **RunQSpinLib_tdvp.m** secondly to calculate the ground-state spin dynamics.
Run **PlotSkomega.m** to plot the DSS $S(k,\omega)$ of 1D Heisenberg Chain with lattice size $L=64$ and bond dimension $D=200$.

## Hamiltonian ##
$$H=\sum_{i=1}^{L-1} S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + S_i^z S_{i+1}^z$$

## Dynamical Spin Structure factor ##

For the isotropic Heisenberg model, we have
$$S(k,\omega) = 3\int_{-\infty}^{\infty}dt \, e^{i\omega t} 
g(k,t), $$
with the correlation function $g(k,t)=\langle e^{iHt}S^z_{-k}e^{-iHt}S_k^z\rangle$ and $S_k^z = \frac{1}{\sqrt{L}}\sum_m^Le^{imk}S_m^z$. 

Note that $g^*(k,t) =\langle e^{-iHt}S^z_{-k}e^{iHt}S_k^z\rangle =g(k,-t)$, we have 
$$S(k,\omega)=3[\int_{0}^{\infty} dt\, e^{i\omega t}g(k,t) + \int_{0}^{\infty} dt\, e^{-i\omega t}g(k,-t)]=6\int_0^{\infty}dt\, {\rm Re}[e^{i\omega t}g(k,t)].$$

### Window function ### 
In practice, we have a maximum evolution time $t_{\rm max}$ in tdvp. Here we introduce a window function $W(t/t_{\rm max})$ which vanishes with $|t/t_{\rm max}|>1$.

Now we have 
$$S(k,\omega)\simeq 6\int_0^{\rm t_{\rm max}} dt\, {\rm Re}[e^{i \omega t} g(k,t)] W(t/t_{\rm max}).$$

Here we provide two type window function, the Parzen window function and the staircase function, demonstrated in **WindowFunction.m**. Notably, the staircase function may lead to sub peaks, which may confuse the researcher. 

### Energy resolution ###
The finite $t_{\rm max}$ leads to a finite energy resolution. To estimate the energy resolution, we consider the full width at half maximum of the Fourier transform of the window function. In practice, the energy resolution of Parzen window function is $\varepsilon\simeq 8/t_{\rm max}$. Please run **Energyresolution.m** to see the detail.


## DMRG preparation controls

The DMRG `ImportMBSolverPara.m` settings now include:

```matlab
Para.DMRGStepMax = 100; % Budget per stage
Para.DMRGSweepTol = 1e-8; % Relative E and change of maximum bond entropy
Para.DMRGUpdateSites = "auto"; % Or 1 / 2 for a fixed mode
Para.DMRGEnvStorage = 'memory'; % Or 'disk'
Para.DMRGInitialMPSFile = ''; % Or the filename of a previous MPS result
Para.DMRGEnvRoot = fullfile(tempdir, 'qspintn_dmrg');
Para.DMRGKeepEnv = true; % false cleans this run's disk cache on exit
Para.saveMPS = true; % Include T and H for the existing MPS save
```

Auto starts with two-site sweeps and can switch to one-site on an energy rise.
The existing `MCrit`, `DK` and local `tol` values are retained. Check
`Rslt.DMRGInfo.Converged` / `StopReason` and bond convergence before TDVP; a sweep
limit or fixed-two-site energy-rise exit is not convergence. The existing
`MPS/` filename and saved `Para`/`Rslt` format remain compatible with the TDVP
loader; filenames do not distinguish update modes. Historical MAT files may lack
these diagnostics. Disk environment files in `Rslt.EnvDir` are not MPS checkpoints.
The complete interface is documented in the
[DMRG README](../../ManyBodySolver/DMRG/README.md).


## Start from a previous MPS file

Set one parameter in `ImportMBSolverPara.m`, or set it after calling the importer:

```matlab
Para.DMRGInitialMPSFile = 'MPS/previous.mat';
Rslt = GetDMRGRslt(Para);
```

No manual `load` or extraction of the state is needed. The file must contain
`Rslt.T` (save the previous result with `saveMPS=true`) or a standalone `T`.
An empty filename selects random initialization. Current model/field parameters
are used to rebuild the Hamiltonian and environments; the old energy and sweep
history are not reused. Match the old state's lattice size, local dimension,
basis and site ordering. The new `MCrit` must cover every input bond. Two-site or
auto mode can increase the ranks up to the new cap; one-site mode cannot expand
the supplied ranks. The loader/solver checks the tensors, restores right
canonical form and normalizes before solving.

Inspect `Rslt.DMRGInfo.InitialStateUsed`, `InitialMPSFile` and `InitialEnergy` to
confirm the source and its energy under the new Hamiltonian. A new solve still
needs its own convergence checks.

For filename warm starts, `RunQSpinLib_DMRG.m` saves the new MPS with a
`-warmstart-<timestamp>` suffix and preserves the input file. Point the TDVP
loader at that newly saved filename. Cold-start filenames are unchanged.

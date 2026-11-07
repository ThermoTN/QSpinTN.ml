# Example\: Triangular lattice ground-state dynamics
We offer here an example of ground-state spin dynamcis of triangular lattice quantum magnet.

## Hamiltonian ##
$$H=\sum_{\langle i,j \rangle} S_i^\alpha J_{i,j}^{\alpha \beta} S_j^\beta + \sum_{\langle \langle i,j \rangle \rangle } J_{2,xy} (S_i^x S_j^x + S_i^y S_j^y) + J_{2,z}S_i^z S_j^z - \sum_{i=1} (h_x S_i^x + h_y S_i^y + h_z S_i^z),$$
where $\langle i,j \rangle$ labels the nearest neighbor (NN) site, $\langle \langle i,j \rangle \rangle$ labels the next-nearest neighbor (NNN) site, and \
$$
J^{\alpha\beta} = 
\begin{pmatrix}
J_{xy}+2 J_{\rm PD} \cos\varphi
&-2J_{\rm PD} \sin\varphi
&-J_{\Gamma} \sin\varphi \\
-2J_{\rm PD} \sin\varphi
&J_{xy}-2J_{\rm PD} \cos\varphi
& J_{\Gamma} \cos\varphi\\
-J_{\Gamma} \sin\varphi
& J_{\Gamma} \cos\varphi
& J_z
\end{pmatrix}, \nonumber
$$

with $\varphi=\\{0,\frac{2\pi}{3},-\frac{2\pi}{3}\\}$ for three different types of NN bonds.

In **PlotSiteLocation.m**, we provide the lattice geometry, with the red, green, and blue dashed line corresponding to $\varphi=\\{0,\frac{2\pi}{3},-\frac{2\pi}{3}\\}$ respectively.


**Note: We have omitted the parameters that were mentioned earlier.**

Please see **Example_QSPinLib_TriARX_1_tanTRG** and **Example_QSpinLib_SpinChain_4_tdvp**.


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

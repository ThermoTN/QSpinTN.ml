# DMRG

`GetDMRGRslt(Para)` builds the MPO and solves for a ground-state MPS using
one-site, two-site or automatic two-site to one-site sweeps. Both stages support
memory and disk environments. Start with
[Example_SpinChain_3_DMRG](../../Example/Example_SpinChain_3_DMRG/README.md).

## Parameters

`MCrit`, `DK`, `tol` and `DMRGStepMax` are required. The other defaults below also
apply to older callers that omit the new fields.

| Parameter | Meaning / default |
| --- | --- |
| `MCrit` | Maximum MPS bond dimension; one-site sweeps cannot increase initial ranks. |
| `DK` | Maximum local Lanczos dimension; randomly initialized two-site runs use 2 for their first sweep, then the requested value. |
| `tol` | Local Lanczos Ritz-energy change / breakdown tolerance, not an eigenvector residual threshold. |
| `DMRGStepMax` | Maximum complete bidirectional sweeps **per stage**; auto can use up to twice this budget. |
| `DMRGSweepTol` | Default `1e-8`; relative total-energy change and absolute change of maximum bond entropy must both pass. |
| `DMRGUpdateSites` | Default `"auto"`; alternatively numeric `1` or `2` for a fixed mode. |
| `DMRGEnvStorage` | Default `'memory'`; `'disk'` writes individual rolling environment blocks. |
| `DMRGEnvRoot` | Disk cache parent; default `fullfile(tempdir,'qspintn_dmrg')`. Ignored in memory mode. |
| `DMRGKeepEnv` | Default `true`; `false` removes this run's disk directory on return or error. Ignored in memory mode. |
| `DMRGInitialMPSFile` | Default `''`; MAT filename containing a previous `Rslt.T` or standalone `T` for a new run. |
| `saveMPS` | Default absent/false; `true` includes `T` and `H` in `Rslt`. It does not itself write a result file. |

```matlab
Para.MCrit = 100;
Para.DK = 20;
Para.tol = 1e-12;
Para.DMRGStepMax = 100;
Para.DMRGSweepTol = 1e-8;
Para.DMRGUpdateSites = "auto";
Para.DMRGEnvStorage = 'memory'; % Set 'disk' to offload environments
Para.DMRGEnvRoot = fullfile(tempdir, 'qspintn_dmrg');
Para.DMRGKeepEnv = true;
Para.saveMPS = false;
Para.DMRGInitialMPSFile = ''; % Or 'MPS/previous.mat'
Rslt = GetDMRGRslt(Para);
```

## Modes and convergence

Auto starts with two-site sweeps. From the second completed sweep, a resolved
energy rise triggers a switch:

```matlab
Enew - Eold > 16*eps(max([1,abs(Enew),abs(Eold)]))
```

Switching takes priority over convergence. It carries the current MPS and its
retained ranks into one-site sweeps, discards the two-site environments, rebuilds
the one-site right environments, and checks the handoff energy. Sweep numbers
remain continuous. If two-site sweeps converge first, no switch occurs. Fixed
mode `2` retains the historical early exit on an energy increase from sweep 4
onward, reported as `StopReason="energy_rise"`, not convergence.

Both stages call `CheckDMRGSweep` after at least two sweeps in that stage:

```matlab
abs(Enew-Eold) < DMRGSweepTol * max(eps,abs(Eold))
abs(max(Snew)-max(Sold)) < DMRGSweepTol
```

The entropy test compares maximum bond entropies, not the largest per-bond change.
Two-site entropy uses the pre-truncation spectrum from `svdT`; one-site updates
do not truncate. The modes share a stopping rule but have different truncation
semantics.

`DMRGLanczos` uses three-term recursion without reorthogonalization and a small
tridiagonal `eig`. It keeps the source two-site Ritz-energy stopping rule
`(Eritz_old-Eritz_new)/max(eps,abs(Eritz_new))*L < tol/10`, including stopping on
a Ritz-energy rise, and handles exact breakdown at zero energy. One-site updates
fuse the current MPO with the left environment before each local solve and cache
the matrix layouts for repeated applications. Only the current fused block is
held, including in disk mode.

## Initialization and outputs

Cold two-site/auto runs use `InitMPS` for `d=2`. For `d>2` they still require
`InitMPS_d`, which is absent from this repository. One-site mode has its own
right-canonical initializer for `d>=2`, with ranks `min(MCrit,d^i,d^(L-i))`.
Validate each higher-spin model against ED before use. The existing two-site
tensor layout requires at least three sites.

```matlab
[Rslt,T] = GetDMRGRslt(Para);
[T,E0,EE,envDir,info] = DMRGMain(Para,H);
[T,E0,EE,envDir,info] = DMRGMain(Para,H,initialState);
```

To start from an earlier DMRG result, provide its filename:

```matlab
Para.DMRGInitialMPSFile = 'MPS/previous.mat';
Rslt = GetDMRGRslt(Para);
```

The file must contain `Rslt.T` (the usual `saveMPS=true` result) or a standalone
variable `T`; `Rslt.T` takes priority when both exist. `''`, `""` or an absent
field keeps random initialization. Relative filenames resolve from MATLAB's
current directory. The input file is read only; old `Para`, `H`, energy,
environments and sweep counters are not reused. Current `Para` defines the new
Hamiltonian and all numerical controls, so coupling constants and field may change.

The MPS must have the same `L`, `d`, local basis and site ordering as the new
model, and every input bond must be no larger than the new `MCrit`. Tensor shapes,
finite entries and bond connectivity are checked; basis/site semantics remain
the caller's responsibility. Missing files, absent states and incompatible
shapes/ranks cause errors rather than a random-start fallback. A right-to-left
QR pass restores the center to site 1 without SVD truncation and normalizes the
state; its previous global amplitude (`lgnorm`) is discarded. This preparation
also applies to `DMRGMain(Para,H,initialState)` for an in-memory MPS.

A larger new `MCrit` permits two-site/auto to grow the old bonds; fixed one-site
mode stays within the supplied ranks. A smaller `MCrit` than an input bond is
rejected rather than silently truncating it. With supplied input both modes use
the requested `DK` from the first sweep. Sweep counters and environment stores
start afresh: this is a warm start, not checkpoint resumption.

`Rslt` contains total natural-unit energy `E0`, per-site magnetization along the
field `M`, spin operators `Op`, final maximum bond entropy `EE`, and `DMRGInfo`.
At zero field the wrapper still assigns `M=0`; this does not measure spontaneous
order. With `saveMPS=true`, it also includes `T` and `H`. The optional second
output returns `T` regardless of `saveMPS`.

`DMRGInfo` records `Mode`, `TwoSiteSweeps`, `OneSiteSweeps`, `Converged`,
`StopReason` (`converged`, `sweep_limit`, or fixed-two-site `energy_rise`),
`SwitchSweep`, `SwitchDeltaE`, `SwitchEnergy`, and `OneSiteInitialEnergy`.
`InitialStateUsed` reports whether an MPS was supplied, `InitialEnergy` evaluates
the prepared state with the current Hamiltonian before the first sweep, and
`GetDMRGRslt` adds `InitialMPSFile` with the resolved filename (empty for a cold run).
`TwoSiteSweepSeconds` / `OneSiteSweepSeconds` and the corresponding
`*LocalIterations` arrays include every complete sweep, including final and
switching sweeps. Timings include environment I/O, but exclude initialization,
observations and final result saving. Full energy/entropy histories remain in
the log rather than in `Rslt`.

## Local observables

`DMRGObs/ObserveSpin(S,T.A)` returns all local expectations with the operator
convention `S(bra,ket)`. `DMRGObs/Correlate(Si,Sj,T.A)` returns `<Si_i Sj_j>` for
`i<j`, including the final site. Both require normalized MPS tensors right
canonical at sites 2:L, as returned by either sweep mode. Complex operators
such as `Sy` and raising/lowering operators use the physical ket index.
Uncomputed diagonal/lower-triangular correlation entries are `NaN` (previously
zero); do not treat them as measured zeros. The legacy `Observe_ezLv` convenience
wrapper still selects spin-1/2 operators; the shared `MPS2V` helper also assumes
`d=2`. The observable regression uses independent basis enumeration for `d=3`.

## Disk environments and MPS files

Each disk run creates a unique `dmrg_<job-or-local>_<token>` directory beneath
`DMRGEnvRoot`; its absolute path is printed and returned as `envDir` and
`Rslt.EnvDir`. Each `env_XXXX.mat` contains `block`, `slot`, and `side` (`L`/`R`).
Writes use a temporary MAT file followed by replacement. Only the current
stage's rolling environments remain after a switch. Disk mode offloads
environments, not the full MPS/MPO or current Krylov workspace.

With `DMRGKeepEnv=false`, the returned path identifies a directory already
removed. Memory mode creates no directory and omits `Rslt.EnvDir`. Retained
environments alone cannot resume a calculation: they lack a synchronized MPS,
sweep position and complete checkpoint. Explicitly save `Para` and `Rslt` with
`saveMPS=true` for downstream TDVP. The three TDVP preparation examples keep their
existing `MPS/` filenames for cold runs. Filename warm starts append a
`-warmstart-<timestamp>` suffix to the newly saved result. Select that new file
explicitly in the TDVP loader. Update modes do not select separate directories.

## Validation and provenance

The core implementation was ported from the 2026-11-07 Nb3X8 DMRG handoff. The
library port uses a configurable portable cache root and retains this
repository's initializer instead of the Nb3X8 `my_InitMPS` hook. The general
`ObserveSpin`/`Correlate` fixes are also included. Material-specific observable
wrappers, production entry points and historical benchmark files are not included. No speedup for this checkout is inferred from those historical tests.

Local validation used MATLAB R2025b Update 5. It covered independent dense ED
for a complex six-site DM chain, both fixed modes and storage backends, actual
auto handoff, warm starts, sweep budgets, bond caps, cache isolation/cleanup,
nonzero-field magnetization and MPS saving/loading. The four DMRG example
configurations were also compared with ED on six-site systems. Local and
two-site observable contractions were checked against independent state-vector
expectations for complex spin-1/2 and spin-1 states. The maximum configured MPS
bond cap in these checks was 100.

Validation scripts and generated logs are retained locally and are not part of
this update. For a new scientific model, independently compare with ED where
tractable, check bond/Krylov/sweep convergence and verify the requested storage
backend. These small-system checks do not establish large-system convergence or
performance.

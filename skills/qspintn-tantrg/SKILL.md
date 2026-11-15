---
name: qspintn-tantrg
description: Use QSpinTN one-site, two-site or RSVD-CBE tanTRG for thermal MPO evolution, memory or disk environments, and controlled sweep benchmarks.
---

# tanTRG

All source paths are relative to the QSpinTN.ml repository root. Read the
[solver contract](../../ManyBodySolver/tanTRG/README.md) for full parameters,
storage lifetime, normalization and diagnostics.

`ManyBodySolver/tanTRG/GettanTRGRslt.m` builds `InitHam`, initializes `SETTN`,
and calls `tanTRGMain`. Follow `Example/Example_SpinChain_2_tanTRG/` for geometry,
interaction map, natural field and logging. Resolve its specific importer;
similarly named importers configure other solvers.

## Modes and RSVD expansion

`TanTRGUpdateSites="auto"` defaults to the previous `Full2`/`beta_switch` scheme.
Numeric `1` / `2` force fixed one-site / two-site; `"cbe"` enables one-site thermal
projector splitting with RSVD pre-expansion. This is imaginary-time evolution,
not real-time TDVP. Fixed one-site cannot grow the initial MPO ranks.

The [shared CBE contract](../../ManyBodySolver/TensorFunction/CBE.md) defines
`CBERatio`, `CBEOversampling`, `CBEPowerIterations`, `CBESvdTol`, and
`CBEExpansionCap`. Retained bonds are capped at `MCrit`; temporary expansion can
exceed it unless explicitly capped. CBE-tanTRG truncates with `TSRGStol`, updates
the logarithmic norm, and includes backward bond evolution. CBE discarded
weight is not the fixed one-site value zero.

## Storage and outputs

`TanTRGEnvStorage='memory'` is the default; `'disk'` offloads rolling environments
for all three methods. `TanTRGEnvRoot` defaults to `/tmp/ygao`; set the actual
server scratch location when needed. Each run creates an isolated directory.
`TanTRGKeepEnv=true` retains it; `false` cleans on exit/error. `Rslt.EnvDir` is
reported for disk runs. The MPO, SETTN and local Krylov workspace stay in memory.
Environment files alone do not contain a resumable MPO checkpoint.

The result contains `beta`, `LnZ`, `En`, `M`, `EE`, `TrunErr`, final `rho` and
`TanTRGInfo.Sweeps`. Beta starts at `2*tau0` and advances by `4*tau` each complete
sweep; `GetNowtau` selects the step. The endpoint may exceed `beta_max`.
`LnZ` and `En` are per site; the initial energy now uses the same scale/division
as later points. `EE`/`TrunErr` now align with beta, with initial zero placeholders;
historical files may have trimmed diagnostics and an unscaled first energy.

Hooks `GettanTRGCorr` and `savetanTRGRslt` are optional. If present, their errors
propagate; check required observables/output explicitly. Environment storage
and full-result saving are separate. Cm is derived downstream by
`ManyBodySolver/ThDQ_func.m`; see [units](../qspintn-units/SKILL.md).

## Fixed-input checks and timing

`TanTRGSweep(rho,H,Para)` performs exactly one bidirectional sweep for explicit
mode `1`, `2` or `"cbe"` and `Para.tau`. Use normalized right-canonical tensors
(center at site 1), with the overall amplitude in `rho.lgnorm`. Copy the same
input for every method. `TanTRGInfo.Sweeps` and the sweep's second output provide
sweep time (including disk I/O), retained/expanded bonds, truncation and
expansion statistics. Initialization and observables are outside that timer.

Compare small systems with independent ED/matrix exponentials; verify memory
and disk agreement. Then vary retained rank, time step and SETTN precision.
Single-sweep speed ratios do not establish equal accuracy, convergence or
large-system performance. The SVD discarded weight does not include missed
RSVD directions or time-evolution projection error; check smaller time steps
even when that weight is tiny. Derivative heat capacity needs its own convergence
check; finite size and truncation errors are distinct.

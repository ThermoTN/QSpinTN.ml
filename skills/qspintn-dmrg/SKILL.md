---
name: qspintn-dmrg
description: Use QSpinTN one-site, two-site or hybrid DMRG for ground states, filename-based MPS warm starts, memory or disk environments and convergence diagnostics.
---

# DMRG

All source paths are relative to the QSpinTN.ml repository root. Read
[ManyBodySolver/DMRG/README.md](../../ManyBodySolver/DMRG/README.md) for the complete
parameter, stopping and storage contract.

## Configure and run

`ManyBodySolver/DMRG/GetDMRGRslt.m` builds the Hamiltonian with `InitHam` and calls
`DMRGMain`. Use `Example/Example_SpinChain_3_DMRG/` for the interaction map,
geometry, local dimension, natural-unit `Field.h`, and solver inputs. Its importer
sets `MCrit`, `DK`, `tol`, `DMRGStepMax`, `DMRGSweepTol`, `DMRGUpdateSites`,
`DMRGEnvStorage`, `DMRGEnvRoot`, `DMRGKeepEnv`, `DMRGInitialMPSFile`, and `saveMPS` explicitly.

- `DMRGUpdateSites="auto"` is the default: start with two-site sweeps and switch
  to one-site after a resolved energy increase between completed sweeps. Numeric
  `1` / `2` select fixed modes. Auto can finish without a switch.
- `DMRGStepMax` is a separate budget per stage. `DMRGSweepTol` defaults to `1e-8`
  and tests relative energy change and change of maximum bond entropy. `tol`
  governs local Lanczos Ritz-energy change/breakdown, not a residual norm.
- `DMRGEnvStorage='memory'` is the default. `'disk'` offloads rolling environments;
  the MPS/MPO and local workspace remain in RAM. `DMRGEnvRoot` defaults to
  `fullfile(tempdir,'qspintn_dmrg')`. Each run creates a unique subdirectory.
  `DMRGKeepEnv=true` retains it; `false` cleans it on exit.

Cold two-site/auto runs initialize `d=2` with `InitMPS`. For `d>2` they still call
the missing `InitMPS_d`; do not launch them without a compatible initializer or
valid supplied state. One-site mode has a generic `d>=2` right-canonical
initializer and ranks capped by `MCrit` and the Hilbert space of each cut. It
cannot increase those ranks. The two-site layout requires `L>=3`.

For a previous DMRG MPS, set `Para.DMRGInitialMPSFile='MPS/previous.mat'` and
call `GetDMRGRslt(Para)` normally. The MAT file must contain `Rslt.T` or standalone
`T` (`Rslt.T` takes priority). Set the filename after the example importer or edit
its default. Empty/absent selects a cold start; relative paths use the MATLAB
working directory. Missing/incompatible input is an error, never a random fallback.
The previous run needs `saveMPS=true` when saving `Rslt`.

Current `Para` defines the new Hamiltonian, fields and numerical controls; no
old Hamiltonian, energy, environments or sweep counters are reused. Match `L`,
`d`, local basis and site ordering. The solver checks shapes/connectivity/finiteness
and requires each input bond to fit `MCrit`, then restores right canonical form
with QR and normalizes without SVD truncation. This also applies to an in-memory
MPS passed to `DMRGMain(Para,H,initialState)`. A larger cap allows two-site/auto
bond growth; fixed one-site sweeps cannot expand the supplied ranks. A smaller
cap than the old bonds is rejected. Supplied states use full `DK` from the first
sweep; cold two-site initialization uses `DK=2` initially. This is a new solve,
not checkpoint resumption.

## Read outputs and prepare dynamics

`[Rslt,T]=GetDMRGRslt(Para)` returns total energy `E0`, per-site field-direction
magnetization `M`, spin operators `Op`, maximum bond entropy `EE`, and `DMRGInfo`.
`saveMPS=true` adds `Rslt.T` and `Rslt.H`; it does not write a file. The optional
second output returns `T` independently of that flag. At zero field the wrapper
assigns `M=0`; this does not establish absence of spontaneous or local order.

Inspect `DMRGInfo.Converged` and `StopReason`, both stage counts and switch/handoff
energies. `InitialStateUsed` confirms whether a supplied state was used,
`InitialMPSFile` records its resolved input filename, and `InitialEnergy` is its
energy under the current Hamiltonian before the first sweep. `sweep_limit` and fixed-two-site `energy_rise` are not convergence.
Per-stage `*SweepSeconds` and `*LocalIterations` include all complete sweeps but
exclude initialization/observations/final saving. Full energy and entropy
histories are printed in the log. `EE` uses a pre-truncation spectrum in two-site
mode and an untruncated spectrum in one-site mode.

Disk runs return `Rslt.EnvDir`; with `DMRGKeepEnv=false` it has already been
removed. These rolling environments contain neither the MPS nor a sweep
checkpoint and cannot alone resume a run. Memory runs omit this field.

For [TDVP](../qspintn-dynamics/SKILL.md), save `Para` and `Rslt` with `saveMPS=true`
and verify `T`, `H`, total `E0`, and diagnostics before reuse. The three DMRG-to-TDVP
examples append `-warmstart-<timestamp>` to newly saved results when an input MPS
file is supplied, preserving the old file. Select the new filename in the TDVP
loader. Filenames still do not encode update mode.
Historical MAT files may lack the new diagnostics.

For local order and correlations, `DMRGObs/ObserveSpin` and `Correlate` use
`S(bra,ket)` and expect a normalized MPS right canonical at sites 2:L.
`Correlate` only evaluates `i<j`; diagonal/lower entries are `NaN`, not measured
zeros. The `Observe_ezLv` convenience wrapper and shared `MPS2V` helper remain
spin-1/2-specific.

## Verification and interpretation

For a small tractable system, compare with independent ED and exercise both
storage backends, mode switching, cache cleanup and MPS saving/loading. Validation
scripts/logs for this update are local development artifacts, not distributed
entry points. Compare each scientific model with ED when tractable, increase
bond dimension/Krylov budget, and check independent seeds for competing states. Sweep convergence and one-site zero
truncation error do not establish exactness. Do not extrapolate another project's
benchmark to this checkout or a different system size. DMRG is a ground-state
solver, not a direct finite-temperature Cm or susceptibility solver; geometry
width, ordering and finite size still matter.

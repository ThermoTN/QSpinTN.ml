---
name: qspintn-tantrg
description: Use QSpinTN tanTRG for imaginary-time thermal MPO evolution and assess step schedules, truncation, and derived heat capacity.
---

# tanTRG

All source paths below are relative to the QSpinTN.ml repository root.

## Entry points and inputs

`ManyBodySolver/tanTRG/GettanTRGRslt.m` runs `InitHam`, `SETTN`, `InitTSRGEnV`, `twosite_tanTRG`, and `onesite_tanTRG`. This is imaginary-time thermal evolution, distinct from real-time TDVP dynamics.

Follow `Example/Example_SpinChain_2_tanTRG/` for model, interaction map, geometry, natural field, and logging setup. Its `Example/Example_SpinChain_2_tanTRG/ImportMBSolverPara.m` provides `tau0`, `beta_max`, `tau_step`, `beta_c`, `TSRGStol`, `Full2`, `beta_switch`, `MCrit`, and SETTN/variational controls. Resolve this specific importer on the MATLAB path; similarly named importers configure different algorithms.

## Outputs and schedules

The raw result includes `beta`, `LnZ`, `En`, `M`, `EE`, `TrunErr`, and the final thermal MPO `rho`. Evolution starts at `beta=2*tau0`; each complete sweep advances beta by `4*tau`, with `ManyBodySolver/tanTRG/GetNowtau.m` selecting the schedule. `Full2` and `beta_switch` control the transition from two-site growth to one-site evolution. Inspect achieved endpoints; the last step can exceed the requested maximum.

`GettanTRGRslt` trims the diagnostic `EE`/`TrunErr` endpoints, so these arrays need not match the thermal grid. Record the actual sampling coordinates. In the physical interface `ManyBodySolver/GetResult.m`, Cm is derived from energy by `ManyBodySolver/ThDQ_func.m` using `csaps` and interior beta coverage; it is not a raw tensor observable on every saved beta. See [units](../qspintn-units/SKILL.md).

## Convergence and hooks

Benchmark small finite geometries against ED, then vary `MCrit`, time-step schedule, SETTN initialization, and the two-site/one-site transition. Check Cm independently because numerical differentiation amplifies energy errors. Finite-size and MPO truncation errors are separate convergence questions.

Evolution calls `GettanTRGCorr` and `savetanTRGRslt` inside a catch block. Select or implement hooks for the desired observables and output location; missing hooks or failed writes can become warnings while evolution continues. Verify required files and observables explicitly. Hook availability and solver documentation do not by themselves provide a downstream execution adapter.

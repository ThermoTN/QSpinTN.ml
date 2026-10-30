---
name: qspintn-iltrg
description: Use QSpinTN iLTRG for repeated two-site spin-chain thermodynamics and assess Trotter, bond-dimension, and output-grid limitations.
---

# iLTRG

All source paths below are relative to the QSpinTN.ml repository root.

## Entry points and inputs

`ManyBodySolver/iLTRG/GetiLTRGRslt.m` calls `ManyBodySolver/iLTRG/RG_func/iTEBD.m` for a repeated two-site tensor network. This is the chain/Trotter implementation; an arbitrary finite-cluster interaction map is not sufficient input.

Use `Example/Example_SpinChain_6_iLTRG/RunQSpinLib_iLTRG.m` as the direct-call interface: `Trotter_Name`, `Model`, `d`, natural-unit `Field.h`, `ThDQ='Cm&Chi'`, `WorkingMode`, `ManyBodySolver`, and `fileID`, then `GetiLTRGRslt(Para, Para.ThDQ)`. Its `Example/Example_SpinChain_6_iLTRG/ImportMBSolverPara.m` supplies `tau`, `N_max`, `TroOrd`, `D_max`, and `InterNum`. Verify which example-local importer MATLAB resolves. The example is itself marked not recommended; treat it as an interface reference requiring validation.

`ManyBodySolver/iLTRG/InitHam_trotter.m` constructs the alternating two-site Hamiltonians. Validate the Trotter terms and field convention before interpreting a different model or anisotropic g tensor.

## Grids, outputs, and convergence

`ManyBodySolver/iLTRG/ThDQ/Tem_cal.m` uses `beta_i = 2*tau*i`. Reducing `tau` without increasing `N_max` loses low-temperature coverage. Compare multiple `tau` and `D_max` values at a common physical temperature window. The available example uses first-order Trotter evolution.

Heat capacity comes from derivatives of `LnZ`; endpoints are removed and `T`, `Cm`, and `M` are interpolated again. `LnZ` is retained on the earlier grid and must not be paired blindly with the final `T`. Inspect array lengths and grids explicitly. Direct outputs are in natural units; plotting `M/norm(Field.h)` as in the example does not itself produce molar SI susceptibility.

## Current implementation limits

- A pure `'Chi'` request leaves `LnZ_l` undefined in `GetiLTRGRslt`; the example's `'Cm&Chi'` route avoids that branch. Do not assume the generic `GetResult(..., 'Chi')` route works for iLTRG.
- In `InitHam_trotter`, the right-operator branches for `'Sz2'` and `'Id'` assign `Sl` instead of `Sr`. Models using these terms require a code correction and an independent Hamiltonian check before numerical use.
- The direct example uses a vector natural field; the shared physical configuration uses field magnitude and direction. Trace the selected route before converting units or g factors.

Check a known chain limit or controlled comparison with increasing finite ED chains, accounting for the different boundary conditions and thermodynamic limit.

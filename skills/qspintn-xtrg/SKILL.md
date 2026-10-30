---
name: qspintn-xtrg
description: Use QSpinTN XTRG for finite-MPO thermodynamics, including SETTN initialization, temperature grids, and convergence checks.
---

# XTRG

All source paths below are relative to the QSpinTN.ml repository root.

## Entry points and inputs

`ManyBodySolver/XTRG/GetXTRGRslt.m` builds the finite Hamiltonian MPO with `InitHam`, calls `ManyBodySolver/XTRG/XTRG_Main.m`, initializes through `ManyBodySolver/TensorFunction/SETTN.m`, and evolves with `ManyBodySolver/XTRG/ExpEvo_conj.m`.

Follow `Example/Example_SpinChain_5_XTRG/` for interaction-map/model/geometry/field setup and its `Example/Example_SpinChain_5_XTRG/ImportMBSolverPara.m`. Solver controls include `tau0`, `It`, `MCrit`, `D_list`, `sweep_step`, `VariProd_step_max`, `VariSum_step_max`, and `SETTN_init_step_max`. `MCrit` controls initial compression while `D_list` sets evolution bond dimensions; they are not interchangeable. Resolve the correct example-local importer on the MATLAB path.

Use `[raw, physical] = GetResult(conf, K_min, ThDQ)` for the package's physical thermal output, with a valid package configuration. `K_min` is the lowest requested temperature in kelvin; `ThDQ` is `'Cm'` or `'Chi'`. Model parameters belong in `conf`, not in the second argument. See [units](../qspintn-units/SKILL.md).

## Grids and outputs

For Cm the entry point interleaves starts `ST=[0,1,2,3]`; for Chi it uses `[0,2]`. Each starts at `tau0*2^(ST/4)` and evolves by inverse-temperature doubling. The combined result sorts `beta`, `LnZ`, `M`, `En`, `EE`, and `TrunErr` together.

Physical Cm is produced by `ManyBodySolver/ThDQ_func.m`, which fits energy against log(beta) using `csaps`, differentiates the spline, and drops four points at each boundary when selecting the output window. The resulting Cm grid differs from the raw energy grid (`betaOri` is retained by `GetResult`). Reserve endpoint coverage and check returned `T_l` before comparison with data. Smooth energy alone does not establish converged heat capacity.

## Applicability and verification

Finite chains, strips, or clusters require a correct interaction map and MPO ordering. The method does not automatically remove finite-size or finite-width effects. Compare against ED on exactly the same small Hamiltonian, then vary evolution bond dimensions, initialization accuracy, and sweeps separately.

Check MATLAB availability of `csaps` for physical Cm. `ManyBodySolver/TensorFunction/svdT.m` normally uses MATLAB `svd` and may fall back to `svd_lapack`; a bundled binary must match the running platform. Archive temperature coverage and truncation diagnostics with the result, and separate numerical error from changes caused by increasing geometry size.

---
name: qspintn-dynamics
description: Use QSpinTN real-time TDVP for ground-state correlation dynamics, including MPS preparation, momentum operators, time units, and saving hooks.
---

# Real-time TDVP dynamics

All source paths below are relative to the QSpinTN.ml repository root.

## Preparation and entry point

`ManyBodySolver/tdvp/GettdvpRslt.m` takes `(Para, rho, H)` and has no returned result argument. Its inputs are a ground-state MPS and Hamiltonian MPO, typically `T` and `H` saved by [DMRG](../qspintn-dmrg/SKILL.md). Supply the corresponding total energy as `Para.E0`; all must describe the same Hamiltonian, field, and geometry.

`Example/Example_SpinChain_4_tdvp/RunQSpinLib_tdvp.m` demonstrates the call. Its `Example/Example_SpinChain_4_tdvp/InitKSpaceIntr.m` constructs `Para.QuIntr` from `DY.k` and `DY.QulocOp`, using the chain operator `Q = sum_j exp(i*k*j)*O_j/sqrt(L)`. Other geometries require their corresponding positions and momentum convention; a chain helper is not universal.

Set `tau`, `tmax`, `MCrit`, `Stol`, `Full2`, `t_switch`, and a valid `fileID`. Distinguish truncation control from the internal Lanczos exponential controls. The chain example loads a specific ten-site MPS under its `MPS/` directory; that file is tracked in this checkout. Verify its saved `Para`, `T`, `H`, and `E0` before reuse. For a changed model, field, geometry, or accuracy target, prepare a matching state with `Example/Example_SpinChain_4_tdvp/RunQSpinLib_DMRG.m` and its importer, which sets `saveMPS=true`, then select the new file explicitly.

## Observable and output contract

The driver evolves `Q|psi0>` and records the complex correlation `C(t) = <psi0|Q^dagger exp(-i*H*t) Q|psi0>*exp(i*E0*t)`. Saved results contain `time_list`, complex `Obs`, `TrunErr`, and `EE`. The output interval is `2*tau`; `tmax` may be exceeded by a final interval. Time is in inverse Hamiltonian energy units with hbar=1, not automatically picoseconds.

Both evolution stages call `savetdvpRslt(Para, Rslt)` within a catch block. Provide a compatible saving hook on the MATLAB path and check its actual output location; examples supply hooks in their own directories. The chain hook writes under `Rslt_DY/` relative to the working directory using `GetFileName(Para)` and does not create that directory. Create the chosen output directory or supply a hook that writes to the intended workspace. Because the driver has no return argument, a missing or failed hook may leave no usable result even when evolution runs. A documented solver is not automatically an executable downstream adapter.

## Verification and interpretation

Check the t=0 equal-time value, normalization and energy conservation, and a small-system ED correlation. Vary time step, bond dimension, total time, and two-site/one-site transition. A one-site `TrunErr=0` only means that stage records no truncation loss; it does not certify zero integration or projection error.

For a frequency spectrum, explicitly choose Fourier sign, normalization, connected/disconnected contribution, window, and broadening. Finite-time resolution and finite-size spectral structure remain physical interpretation limits, even after time evolution converges. This entry point describes ground-state real-time correlations, not a general finite-temperature dynamical solver.

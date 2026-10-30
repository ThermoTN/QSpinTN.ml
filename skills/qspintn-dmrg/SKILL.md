---
name: qspintn-dmrg
description: Use QSpinTN DMRG for ground-state energy, magnetization, and saved MPS preparation, with current initialization and output limitations.
---

# DMRG

All source paths below are relative to the QSpinTN.ml repository root.

## Entry points and inputs

`ManyBodySolver/DMRG/GetDMRGRslt.m` constructs the Hamiltonian with `InitHam` and calls `ManyBodySolver/DMRG/DMRGMain.m`. Follow `Example/Example_SpinChain_3_DMRG/` for `IntrcMap_Name`, `Model`, `d`, geometry, `L`, natural-unit `Field.h`, `ManyBodySolver`, and `fileID`. The local `Example/Example_SpinChain_3_DMRG/ImportMBSolverPara.m` supplies `MCrit`, Krylov dimension `DK`, local tolerance `tol`, and `DMRGStepMax`.

The current `DMRGMain` initializes `d=2` through `ManyBodySolver/DMRG/Tensor_Func/InitMPS.m`, which uses random states. For `d!=2` it calls `InitMPS_d`, but that function is absent from this checkout. A higher-spin run requires a compatible initializer and validation before execution; there is no current `my_InitMPS` fallback.

## Outputs and physical interpretation

`GetDMRGRslt` returns `E0`, `M`, and local spin operators `Op`. Set `saveMPS=true` to also return `T` (MPS) and `H` (MPO), needed by [dynamics](../qspintn-dynamics/SKILL.md). `E0` is the total natural-unit energy, while `M` is divided by `L`. The wrapper assigns `M=0` when the applied field is zero; this does not prove absence of spontaneous or local order.

The present wrapper does not return `EE`, sweep history, or a convergence status. Although `DMRGMain` prints energy, entropy, and truncation information, those diagnostics must be captured or exposed separately if needed.

## Convergence and limits

Check ground-state energy against ED for a matching small Hamiltonian and test increasing bond dimension, Krylov dimension, and sweeps. Record the random seed and compare independently initialized states when competing solutions are plausible. The global sweep termination uses hard-coded energy/entropy thresholds and can also stop after an energy increase; `Para.tol` governs the local Krylov calculation, not all stopping conditions.

DMRG computes a ground state, not finite-temperature Cm or susceptibility. A two-dimensional interaction map still yields a finite MPO geometry whose width and ordering affect cost and accuracy. Save `Para`, `T`, `H`, and total `E0` together if preparing a TDVP run; the field, geometry, and energy convention must match.

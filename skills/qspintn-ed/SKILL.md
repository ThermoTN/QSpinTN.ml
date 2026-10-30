---
name: qspintn-ed
description: Use QSpinTN exact diagonalization for small finite spin systems, thermal benchmarks, and high-temperature model exploration.
---

# Exact diagonalization

All source paths below are relative to the QSpinTN.ml repository root.

## Entry points and inputs

Use `ManyBodySolver/ED/GetEDRslt.m` for the ordinary ED implementation. It builds the Hamiltonian through `GetLocalSpace`, `AutomataInit`, `AutomataInit1Site`, and `MPO2H`. `ManyBodySolver/ED/ED_func/ED_Cm.m` diagonalizes the full matrix; the susceptibility branch also needs eigenvectors. Estimate memory using Hilbert dimension `(2S+1)^L`, including dense matrices and workspace.

Follow `Example/Example_SpinChain_1_ED/` for the direct `Para` interface: interaction map, model parameters, local dimension `d`, geometry and boundary conditions, `L`, field, `beta_list`, and `ThDQ='Cm'` or `'Chi'`. Example-local `ImportMBSolverPara` functions share a name; check MATLAB's `which` to ensure the intended one resolves.

For experimental units, `ManyBodySolver/GetResult.m` exposes `[raw, physical] = GetResult(QMagenConf, K_min, ThDQ)`, where `K_min` is the lowest requested temperature in kelvin. `Example/Example_QMagen_SpinChain/RunOpt_ED.m` illustrates the configuration structure; `UtilityFunc/PassPara.m` and `UtilityFunc/GetModel.m` define how it reaches the solver. The direct ED example uses a different `Para` interface and does not apply the physical conversion. See [units](../qspintn-units/SKILL.md).

For a direct experimental calculation, set `Para.Model.ES` in kelvin and express the interaction-map couplings in units of that scale. Add `UtilityFunc/` to the MATLAB path, then specify the physical field and moment operator explicitly, for example for isotropic g:

```matlab
Para.Field.B = [0, 0, B_T]; % B_T is nonzero, in tesla
Para.UnitCon = GetUnitCon(Para);
Para = GetField(Para); % scalar Field.h and unit-vector Field.hdir
Para.g_fec.g = [g, g, g];
Para.g_fec.dir = {[1,0,0], [0,1,0], [0,0,1]};
```

The direct example omits `g_fec`, so `AutomataInit1Site` defaults to g=1 and accepts its vector natural field; adding `g_fec` also requires the scalar `Field.h` and `Field.hdir` convention above. For three spin-1/2 sites set both `Para.L` and `Para.Geo.L` to 3, `d=2`, and choose the intended boundary condition. Run `ThDQ='Cm'` and `'Chi'` separately. On the direct result's own grid, `T=ES./raw.beta`, `Cm=raw.Cm*Para.UnitCon.Cm_con`, and `chi_CGS=(raw.M/Para.Field.h)*Para.UnitCon.Chi_con/(4*pi)`. This last expression applies only to the Chi result with a nonzero field and the moment convention above; convert experimental uncertainties by the same unit factor as their values.

## Outputs and numerical interpretation

- Cm: `H`, `Op`, `beta`, `T`, `Cm`, `En`, `F`; heat capacity, energy, and free energy are per site.
- Chi: `H`, `Op`, `Mop`, `beta`, `T`, `M`; `M` is per site and the field-dependent moment operator includes the model's g factor.
- Natural inverse temperature is `beta`; the physical temperature is `ES/beta` kelvin.

ED is exact for the specified finite Hamiltonian at low as well as high temperature. Using a small cluster as a material approximation requires finite-size checks, particularly at low temperature. Compare spectra and thermal observables on the same geometry before using ED to benchmark a tensor method.

`GetResult` resamples susceptibility onto its own beta grid with a spline. Inspect the returned `T_l` coverage rather than assuming `beta_list` determines every physical output point; do not extrapolate across missing coverage. Its Cm conversion catches errors and can substitute zeros, so a zero curve alone does not establish a valid calculation.

`ED_C`, `ED_SM`, and `ManyBodySolver/ED/ED_func/ED_M.m` are different interfaces, not drop-in thermal replacements. In particular `ManyBodySolver/ED/ED_func/ED_M.m` writes a benchmark file and enters `keyboard`; do not use it as an unattended production entry point.

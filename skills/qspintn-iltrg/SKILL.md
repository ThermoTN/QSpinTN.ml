---
name: qspintn-iltrg
description: Use QSpinTN iLTRG for repeated two-site spin-chain thermodynamics and assess Trotter, bond-dimension, and output-grid limitations.
---

# iLTRG

Source paths are relative to the **QSpinTN.ml repository root**.
For standalone thermal calculations, prefer `ManyBodySolver/iLTRG/iLTRGForward.m`:
`[result,state] = iLTRGForward(request,options)`. The numerical API version is
`qspintn.iltrg.forward.v1`; check
[interface.json](../../ManyBodySolver/iLTRG/interface.json) before integrating it.
Read [the API guide](../../ManyBodySolver/iLTRG/README.md) for fields, units,
callbacks, direct MATLAB examples and checkpoint handoff.

## Applicability and model preparation

The physical forward profile is validated for spin-1/2, an infinite chain with a
repeated two-site cell, and real nearest-neighbor bilinear exchange terms.
An arbitrary finite interaction map or a two-dimensional lattice is not input
for this profile. Choose another method when these assumptions do not apply.

Supply `terms_ab` and `terms_ba` as operator pairs (`x/y/z`) with `coefficient_K`.
BA defaults to AB. Both bonds are ordered left-to-right; uniform Dz uses `xy=+Dz`
and `yx=-Dz` on both. Keep `field_T` and diagonal `g` separate from exchange terms.
Validate a new Hamiltonian against a small independent matrix before fitting.
`Chi` requires a single nonzero principal-axis field and means finite-field M/H.
The mere availability of this solver does not enable a consuming application's
fitting adapter; check its supported execution tools separately.

## Numerical choices

- `tau` is in K^-1 and `beta_n=2*n*tau`; `energy_scale_K` is handled internally.
  When reducing tau, increase N_max to retain the low-temperature coverage.
- Cm uses a centered lnZ derivative. Keep an extra step beyond the lowest
  desired T. Missing derivative endpoints are omitted, never extrapolated.
- Compare multiple D_max and tau values at common physical temperatures.
  Varying both together measures combined change, not separate error sources.
- Inspect transfer residuals, actual bond dimensions, cutoff warnings and
  negative outputs. Do not silently smooth or clip numerical defects.
- This implementation uses the existing first-order Trotter sequence.
  Small-D unit tests are not proof of convergence in the target low-T regime.

## Outputs and continuation

`result.data` has explicit experimental units: Cm in J/(mol K), Chi in
emu/mol/Oe, per mol spin. `result.raw` also keeps natural observables, μB/spin
moments, and explicitly named `Chi_qspintn_native`. The native convention is
4π times CGS; the existing thermal-unit contract is unchanged. Read the
[units skill](../qspintn-units/SKILL.md) when combining these with other sources.
Conversions belong to `UtilityFunc/ConvertThermalObservables.m`; shared constants
are in `QSpinThermalUnits.m` and preserve `GetUnitCon` values.

Use `options.checkpoint_path` for file restart or `options.state` for an in-memory
state. `pause_after_steps` or `should_stop` can pause between complete iterations;
`on_progress` reports state to the caller. The application chooses paths, budgets,
threads and whether to resume; the solver verifies model/settings/source identity.
Changed numerical conditions require a fresh state. The API never creates
application directories or calls an LLM.

## Natural-unit interfaces

Existing `GetiLTRGRslt`/`iTEBD` now share `iLTRGEvolve`; observation wrappers share
`iLTRGObserve`. Pure Chi works, all result arrays use their native common grid,
and no lower-temperature interpolation is added. `InterNum` no longer controls
this wrapper's grid. Bare spin projections and natural temperatures are not the
physical API's experimental units. The right-side Id/Sz2 assignment in
`InitHam_trotter` is corrected; that does not expand the validated physical
forward profile to arbitrary spin or interaction classes.

Numerical validation lives in `tests/test_iltrg_forward.m`: free spins, independent
DM dimers, unit/energy scaling, observable normalization and exact restart.

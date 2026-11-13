---
name: qspintn-iltrg
description: Use QSpinTN iLTRG for spin-1/2 nearest- and next-nearest-neighbor chain thermodynamics, including J1-J2 and directed DM terms, and assess grouping, Trotter, bond-dimension and output-grid limitations.
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

The physical forward profiles use spin-1/2, an infinite chain with a repeated
physical A/B cell, and real bilinear exchange through nearest and uniform
next-nearest neighbors. The range-two profile groups two physical spins into
each local tensor; the physical spin is still 1/2.
An arbitrary finite interaction map or a two-dimensional lattice is not input
for this profile. Choose another method when these assumptions do not apply.

Supply `terms_ab` and `terms_ba` as operator pairs (`x/y/z`) with `coefficient_K`.
BA defaults to AB. Both bonds are ordered left-to-right; uniform Dz uses `xy=+Dz`
and `yx=-Dz` on both. Keep `field_T` and diagonal `g` separate from exchange terms.
Validate a new Hamiltonian against a small independent matrix before fitting.
`Chi` requires a single nonzero principal-axis field and means finite-field M/H.
The mere availability of this solver does not enable a consuming application's
fitting adapter; check its supported execution tools separately.

For J1-J2, put the nearest-neighbor tensor in `terms_ab`/`terms_ba`, and the
uniform distance-two tensor in optional `terms_nnn`. All three lists use the
same operator-pair/coefficient_K structure. Isotropic J2 requires xx, yy and zz
with the same coefficient. `terms_nnn` acts on both A-A and B-B bonds i→i+2;
directed DM keeps its sign on each physical left-to-right bond. Separate J2
tensors on the two sublattices and bonds beyond distance two are unsupported.
Use exactly the `operators` and `coefficient_K` fields per term. Geometry comes
from the enclosing list; extra fields such as `distance`, `sites` or
`sublattice` are rejected and cannot extend that geometry.

Empty/omitted `terms_nnn` uses the original local-dimension-2 path. A nonempty
list, including explicitly zero terms, selects local dimension 4 for `(A,B)`.
`terms_ab` is internal to each block; `terms_ba` joins B of one block to A of
the next. Half of each block-local Hamiltonian and field goes on each adjacent
block bond; the NNN tensor contributes both A-A and B-B across blocks. Do not
add field energy again as an exchange term. The evolution cell then contains
four physical spins. This construction is documented in the API guide and
does not convert arbitrary finite interaction maps into infinite-chain input.

The `forward.v1` extension is additive: inspect `interface.json` profiles for
range-two capability. autoQMagen's existing formal fitting adapter is still
nearest-neighbor-only until its own ModelSpec mapping and capability checks are
updated and tested. A skill or a newer QSpinTN checkout alone does not enable it.

## Numerical choices

- `iLTRGForward(request)` takes physical `request.tau` in K^-1:
  `beta_n=2*n*request.tau`. Its internal `Para.tau` is dimensionless and equals
  `request.tau*energy_scale_K`. If a caller exposes dimensionless tau, convert it
  once at this API boundary: `request.tau=tau_dimensionless/energy_scale_K`.
  The direct API accepts an explicit positive energy_scale_K for each request;
  it does not choose an optimizer policy or shrink request.tau automatically.
  A fitting application may derive E0 per candidate from the largest absolute
  physical interaction coefficient and largest Zeeman field energy across its
  evaluated curves, using one E0 for that candidate. Expand ratios/factors first.
  Keep E0 and the resolved physical step unchanged within a tensor continuation.
  For example, actual dimensionless tau=0.01 and E0=2 K mean request.tau=0.005 K^-1.
  Do not multiply or divide twice. Reducing the physical step requires increasing
  N_max to retain low-T coverage; changing E0 while holding physical tau and
  physical coefficients fixed is only a change of internal energy representation.
  Include physical J2 coefficients when a caller computes E0. Grouping introduces
  no extra factor of two in E0, beta or request.tau.
- Cm uses a centered lnZ derivative. Keep an extra step beyond the lowest
  desired T. Missing derivative endpoints are omitted, never extrapolated.
  Its highest finite native point is `1/(4*request.tau)` K, equivalently
  `energy_scale_K/(4*tau_dimensionless)`. Raw finite-field Chi starts one step
  earlier, at `1/(2*request.tau)` K; a consuming fitter may enforce a common
  narrower Cm-compatible window, so check its contract rather than assuming it.
- autoQMagen fitting treats dimensionless tau as an approved upper bound (default
  0.01, D=50). If the candidate grid cannot cover Tmax, its adapter uses
  tau_effective=coverage_safety_factor*E0/(4*Tmax), default factor 0.99, then passes
  tau_effective/E0 to this API. Otherwise it keeps the requested bound. The
  application checks low-T steps/budgets before evolving all curves, records
  actual values, and does not crop data. These are consuming-application rules,
  not behavior of a standalone iLTRGForward call.
  Its coverage preflight follows the API's actual multiply/divide round trip;
  at a floating-point endpoint it can reduce the input step by a few ULPs,
  recording that adjustment instead of allowing extrapolation.
- Compare multiple D_max and tau values at common physical temperatures.
  Varying both together measures combined change, not separate error sources.
  The d=4 block path costs more and truncates across different cuts; do not
  transfer an NN convergence claim solely because D_max has the same value.
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

Grouped lnZ and Cm are already per physical spin; the bare moment uses
`(S_A+S_B)/2` and stays bounded by 1/2. Never use a spin-3/2 operator for the
four-dimensional block or divide returned observables by two a second time.

Use `options.checkpoint_path` for file restart or `options.state` for an in-memory
state. `pause_after_steps` or `should_stop` can pause between complete iterations;
`on_progress` reports state to the caller. The application chooses paths, budgets,
threads and whether to resume; the solver verifies model/settings/source identity.
Changed numerical conditions require a fresh state. The API never creates
application directories or calls an LLM.
Changing J2 or switching the grouping path invalidates a tensor checkpoint.

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
`tests/test_iltrg_j1j2.m` adds independent physical-ring matrix construction,
grouped free spins, exact classical J1z–J2z transfer results, high-T quantum
J1–J2 comparisons against two ED sizes, separate D/tau comparisons, J1=0 and
J2=0 limits, scale invariance and grouped checkpoint checks. Keep artifacts
outside the source checkout. These tests do not establish low-T convergence.

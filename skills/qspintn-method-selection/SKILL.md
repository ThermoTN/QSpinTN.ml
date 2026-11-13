---
name: qspintn-method-selection
description: Select an appropriate QSpinTN many-body method from observable, geometry, temperature, and convergence needs.
---

# Select a many-body method

All source paths in this skill collection are relative to the QSpinTN.ml repository root. Method descriptions expose package knowledge; an application must separately implement its execution adapter.

| Task | Candidate | Principal check |
| --- | --- | --- |
| Small finite Hamiltonian, thermal curves, reference spectrum | [ED](../qspintn-ed/SKILL.md) | Full Hilbert-space memory and finite-size effects |
| Spin-1/2 chain with NN and uniform NNN bilinear exchange | [iLTRG](../qspintn-iltrg/SKILL.md) | Physical two-spin grouping for J1-J2, field/unit normalization, time-step and bond convergence |
| Finite chain, strip, or cluster thermal MPO | [XTRG](../qspintn-xtrg/SKILL.md) | SETTN accuracy, growing bond dimension, derivative/grid coverage |
| Thermal MPO with flexible imaginary-time schedule | [tanTRG](../qspintn-tantrg/SKILL.md) | Step schedule, bond growth, one-site transition, derivative coverage |
| Ground state or MPS preparation | [DMRG](../qspintn-dmrg/SKILL.md) | One/two-site mode, initializer availability, competing states, sweep/bond convergence |
| Ground-state real-time correlations | [TDVP](../qspintn-dynamics/SKILL.md) | Matching ground state, observable/save hook, time and bond convergence |

For material thermodynamics, start with high-temperature data and a tractable finite ED model to constrain energy scales and g factors, then compare larger sizes and suitable thermal tensor methods at lower temperature. Preserve the Hamiltonian, field direction, and [unit convention](../qspintn-units/SKILL.md) across stages. A good high-temperature fit alone does not uniquely determine the effective Hamiltonian.

When comparing solvers, preserve physical K/T inputs while applying the chosen
energy-normalization policy consistently. autoQMagen derives E0 per candidate
from physical couplings and all evaluated field energies; changes between
candidates under the approved rule do not change the model or require a new
strategy. This scale policy does not provide a common time-step definition:
ED has no Trotter tau, iLTRG uses its own centered-Cm coverage condition, and
XTRG/tanTRG retain their native schedules. Read the [units skill](../qspintn-units/SKILL.md)
for physical/native beta and tau boundaries; independent convergence is still needed.

Select based on the requested observable, geometry/interaction map, local spin dimension, field and temperature coverage, numerical precision, and actual runtime resources. Increasing accuracy and increasing system size solve different problems; assess both. On overlapping accessible systems, compare methods before using larger calculations to discriminate model parameters.

The direct iLTRG API accepts uniform next-nearest terms through `terms_nnn` by
grouping two spin-1/2 sites. This is a specific repeated-chain construction,
not generic long-range support. Check the consuming application's adapter:
autoQMagen's existing nearest-neighbor ModelSpec gate is not expanded merely
by updating this package. An explicit finite cluster is still a finite problem,
even when its bonds happen to be nearest or next-nearest neighbors.

Do not select DMRG as a direct finite-temperature Cm solver or real-time TDVP as imaginary-time thermal evolution. Do not infer thermodynamic-limit validity solely from use of a tensor network. Current iLTRG restrictions and DMRG initialization limits are documented in their method skills. DMRG one-site mode has a generic initializer; cold two-site/auto runs still need the absent `InitMPS_d` for higher spin. Memory/disk environments are storage choices, not different physical approximations.

For native optimizer controls and candidate fitting strategies, see [Bayesian optimization](../qspintn-bayesian-optimization/SKILL.md). The consuming application still owns the executed optimizer configuration, residual weighting, data splits, workflow resumption and agent orchestration; method guidance does not change those policies or its supported tools.

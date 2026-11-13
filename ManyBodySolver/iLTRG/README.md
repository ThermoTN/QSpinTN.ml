# iLTRG forward interface

`[result, state] = iLTRGForward(request, options)` performs an infinite-chain
thermal calculation with a repeated physical A/B cell. The numerical API is
`qspintn.iltrg.forward.v1`, declared in [interface.json](interface.json).
Its physical profiles are spin-1/2 with real nearest-neighbor bilinear
coefficients and optional uniform next-nearest-neighbor coefficients, including
directed DM terms. Next-nearest bonds use two-spin blocks internally. This does not implement
arbitrary finite lattices or a fitting optimizer.

## Direct MATLAB use

Add only `ManyBodySolver/iLTRG` to the MATLAB path. The entry point temporarily
loads its own explicit numerical profile and restores the caller's path. No
agent, JSON protocol, process launcher or project-directory layout is needed.

```matlab
addpath('/absolute/path/QSpinTN.ml/ManyBodySolver/iLTRG');
terms = struct('operators', {{'x','x'}}, 'coefficient_K', 1.5);
terms(2) = struct('operators', {{'y','y'}}, 'coefficient_K', 1.5);
terms(3) = struct('operators', {{'z','z'}}, 'coefficient_K', 1.5);
request = struct('spin', .5, 'energy_scale_K', 1, ...
    'terms_ab', terms, 'field_T', [0,0,0], 'g', [1.5,1.5,2], ...
    'observable', 'Cm', 'tau', .01, 'D_max', 32, ...
    'N_max', 53, 'T_min', 1);
result = iLTRGForward(request);
semilogx(result.data.T_K, result.data.values);
ylabel(result.data.unit); xlabel('T (K)');
```

This short example demonstrates the API, not low-temperature convergence.
`terms_ba` defaults to `terms_ab`. AB and BA terms both follow the physical
left-to-right orientation. For uniform `Dz*(Sx_i*Sy_j-Sy_i*Sx_j)`, put `xy=+Dz`
and `yx=-Dz` on **both** bond types. Magnetic fields are specified separately
through `field_T` (T) and diagonal `g`; do not insert the same field twice.
Each term accepts exactly `operators` and `coefficient_K`. The enclosing list
determines its bond: extra fields such as `distance`, `sites` or `sublattice`
are rejected rather than changing the interaction geometry.

## J1–J2 and next-nearest-neighbor terms

Add `terms_nnn`, a struct array with the same `operators` and `coefficient_K`
fields, for bonds from physical site `i` to `i+2`. The same list acts on A–A
and B–B next-nearest bonds. For an isotropic chain,

\[
H/k_B=\sum_i\bigl[J_1\mathbf S_i\cdot\mathbf S_{i+1}
 +J_2\mathbf S_i\cdot\mathbf S_{i+2}\bigr]
 -\frac{\mu_B}{k_B}\sum_{i,\alpha}g_\alpha B_\alpha S_i^\alpha.
\]

Starting from the request above, the following adds isotropic J2=0.4 K:

```matlab
request.terms_nnn = struct('operators', {{'x','x'}}, 'coefficient_K', .4);
request.terms_nnn(2) = struct('operators', {{'y','y'}}, 'coefficient_K', .4);
request.terms_nnn(3) = struct('operators', {{'z','z'}}, 'coefficient_K', .4);
request.D_max = 50;  % starting setting; verify D and tau for the target regime
result = iLTRGForward(request);
```

J1 and J2 may be anisotropic real bilinear tensors; they are not limited to
isotropic exchange. A distance-two DM interaction uses the same directed
`xy=+Dz`, `yx=-Dz` convention in `terms_nnn`. Do not also put a distance-two
term into `terms_ab` or `terms_ba`. The A–A/B–B next-nearest tensors are uniform;
different J2 on these two sublattices, third-neighbor bonds, multi-spin terms,
higher physical spin, and arbitrary non-field on-site terms are outside this API.

Omitting `terms_nnn` or supplying an empty list retains the nearest-neighbor
path with local tensor dimension 2. A **nonempty** list selects the block path,
even if all its coefficients are zero. This explicit zero-J2 path is useful for
checking the grouping convention; it is not expected to have identical finite-D
or finite-step errors to the ungrouped calculation.

Each block contains two physical spin-1/2 sites `(A_k,B_k)` and has local tensor
dimension 4. With dimensionless physical-spin operators `A` and `B`, define
the block-local Hamiltonian `H0` from the `terms_ab` exchange and the field on
both spins. For two adjacent blocks, the constructor uses

\[
\begin{aligned}
H_0 &= \sum_{\alpha\beta}J^{AB}_{\alpha\beta}S_A^\alpha S_B^\beta
      -\mathbf h\cdot(\mathbf S_A+\mathbf S_B),\\
V &= \sum_{\alpha\beta}J^{BA}_{\alpha\beta}S_B^\alpha\otimes S_A^\beta
 +\sum_{\alpha\beta}J^{(2)}_{\alpha\beta}
 (S_A^\alpha\otimes S_A^\beta+S_B^\alpha\otimes S_B^\beta),\\
H_{\rm block,bond} &= V+\tfrac12(H_0\otimes I_4+I_4\otimes H_0).
\end{aligned}
\]

Here all J and h in the numerical constructor are divided by `energy_scale_K`.
Both alternating evolution gates use this same block-bond Hamiltonian; the
transfer cell therefore contains two blocks, or **four physical spins**.
Every physical NN/NNN bond and field appears once in the infinite Hamiltonian.
The block's spin observable is `(S_A+S_B)/2`, not a spin-3/2 `SpinOp(4)`.
`lnZ`, heat capacity and magnetization are normalized per **physical spin**;
do not divide the returned values by two again. `request.spin` remains `0.5`.

The optional field is an additive extension of `forward.v1`; consumers should
inspect the declared capability/profile rather than infer range-two support
from the API version alone. This package extension does not update autoQMagen's
currently nearest-neighbor-only fitting adapter or automatically enable a new
model in another application.

## Temperature and energy scale

All coefficients and `energy_scale_K` are in K. `tau` is a physical inverse
temperature step in K^-1, so `beta_n=2*n*tau` and `T_K=1/beta_n`. Internally
exchange and field energies are divided by the energy scale and the evolution
step is multiplied by it. `N_max` must cover the target temperature and leave
one extra neighbor for the centered heat-capacity derivative.

The two tau conventions are related by
`tau_dimensionless = request.tau * energy_scale_K`.
The low-level `iLTRGEvolve` field `Para.tau` is dimensionless; the public physical
API field `request.tau` remains in K^-1. A caller that exposes dimensionless tau
must supply `request.tau=tau_dimensionless/energy_scale_K`, exactly once.
For example, actual dimensionless tau=0.01 and a reference energy of 2 K mean
`request.tau=0.005 K^-1`. That gives `beta_n=2*n*0.005 K^-1`, and the highest
finite centered Cm point is 50 K. With a 1 K reference it is 25 K.
This direct API accepts the caller's positive `energy_scale_K`; it does not select
an optimizer-wide or candidate-dependent scaling policy by itself. Preserve the
request's scale and physical step when resuming its tensor checkpoint.

The autoQMagen fitting adapter derives E0 separately for each candidate from the
maximum absolute physical interaction coefficient and the largest Zeeman field
energy across all curves in that evaluation. Exactly all-zero interactions and
fields use 1 K; no separate reference is imposed in ModelSpec. The application's
dimensionless tau setting is an upper bound. When E0/(4*tau) does not cover Tmax,
it uses `tau_effective=coverage_safety_factor*E0/(4*Tmax)` (default factor 0.99),
then passes `request.tau=tau_effective/E0`. Otherwise it retains the requested
upper bound. It checks the revised low-temperature step count before numerical
execution and records the actual E0 and step. A standalone call here does not
apply that fitting policy or automatically adapt its physical tau.

Grouping two spins changes the local Hilbert-space size, not the temperature
definition: `beta_n=2*n*request.tau` also holds for J1–J2. There is no extra factor
of two in beta, field, J2 or E0. For a caller-selected energy scale based on the
largest physical coupling, include the distance-two coefficients alongside the
nearest-neighbor ones before normalizing all terms consistently.

In general `Tmax_Cm=1/(4*request.tau)=energy_scale_K/(4*tau_dimensionless)`.
Finite-field Chi has a native point at the first evolution step, hence
`Tmax_Chi=1/(2*request.tau)`. An application may use the shared Cm coverage bound
for joint fitting. These are grid endpoints, not Trotter-error or convergence
guarantees. Changing the reference energy while fixing dimensionless tau changes
the physical step; preserving the same physical calculation requires preserving
`request.tau` and the physical coefficients, as tested by the energy-scale checks.

Optional scientific settings are `transfer_tolerance` (default 1e-12) and
`transfer_max_iterations` (10000). Failed fixed-point convergence, inconsistent
left/right eigenvalues, non-Hermitian input or unphysical spin projection cause
an error, rather than a successful result with an unreliable observable.

## Results and units

| Field | Meaning |
| --- | --- |
| `result.status` | `completed` or `paused`; numerical failures throw an error |
| `result.data` | Ascending `T_K`, `values`, explicit `column` and `unit` |
| `result.raw` | Native beta grid, lnZ per spin, spin projection, natural Cm, converted quantities and transfer diagnostics |
| `result.summary` | API version, model profile, grouping, iteration, resume position, actual bond dimensions, elapsed time and numerical identity |
| `result.hamiltonian` | Actual `Hab`, `Hba` and numerical `Para`; block-path matrices are 16×16 |
| `result.science` / `result.units` | Request, source hashes and physical conventions |
| second return `state` | Tensor evolution state suitable for `options.state` |

`summary.model_profile` is `nearest_neighbor_spin_half` or
`nn_nnn_spin_half_two_site_blocks`. `maximum_neighbor_distance`,
`tensor_local_dimension`, `physical_sites_per_tensor` and
`physical_sites_per_unit_cell` expose the declared profile and grouping; the
last field counts physical spins in the two-tensor evolution/transfer cell.

`UtilityFunc/QSpinThermalUnits.m` centralizes the existing package constants.
`ConvertThermalObservables(raw,field_T,g)` converts per-spin natural quantities:

- `Cm_J_per_mol_K`: J/(mol K), per mole of magnetic spins.
- `M_muB_per_spin`: principal-axis magnetic moment, including g once.
- `Chi_emu_per_mol_Oe`: finite nonzero-field M/H in CGS experimental units.
- `Chi_qspintn_native`: explicitly separate existing SI-cm³/mol convention,
  numerically `4*pi*Chi_emu_per_mol_Oe`.

The existing [thermal-unit contract](../../skills/contracts/thermal_units.json)
is unchanged. Never store the new CGS values under an old native-chi label.
A `Chi` request requires one nonzero principal-axis field. `Cm` also supports
finite fields; for oblique fields the scalar moment/chi conversion is undefined
and returned as NaN, while heat capacity remains available. Zero-field M/H is
undefined, not a zero-field differential susceptibility.

Cm is a centered second difference of lnZ. Raw derivative endpoints are NaN;
`result.data` omits them. No smoothing, extrapolation or clipping of negative
values is performed. This API does not evaluate arbitrary requested temperature
points; consumers must respect the native grid and assess interpolation error.

## Checkpoints and runtime options

The numerical package owns tensor states, compatible restart and atomic `.mat`
checkpoint writes. Callers choose paths, directories, scheduling and resources.
By default there is no checkpoint or artifact output; diagnostics print to stdout.

```matlab
options = struct('checkpoint_path', '/absolute/existing/directory/state.mat', ...
    'checkpoint_interval', 10, 'pause_after_steps', 12);
[partial,state] = iLTRGForward(request, options);
options = rmfield(options, 'pause_after_steps');
continued = iLTRGForward(request, options);  % automatically loads state.mat
% Alternatively: iLTRGForward(request, struct('state', state))
```

The checkpoint parent directory must already exist. Runtime options:

| Option | Purpose |
| --- | --- |
| `checkpoint_path` | Optional absolute state-file path |
| `checkpoint_interval` | Save every N complete iterations, default 10 |
| `pause_after_steps` | Maximum additional steps this invocation, default Inf |
| `state` | In-memory restart; cannot compete with an existing checkpoint file |
| `log_fid` | Caller-owned open MATLAB file handle, default stdout (1) |
| `on_progress` | Function handle receiving a progress struct |
| `should_stop` | No-argument predicate checked between complete iterations |

Final and paused states are saved. Failures retain the last durable complete
checkpoint; partly evolved tensors are not saved as a valid step. State identity
includes scientific settings and numerical source hashes using relative paths.
Changing model, solver settings (including N_max), units or numerical code requires
a new state. Changing only the checkout location does not change source identity.
Thread limits, CPU affinity, process termination, JSON and CSV belong to callers.
Next-nearest coefficients and the physical-spin grouping enter that identity;
a nearest-neighbor tensor state is not a J1–J2 warm start.

## Shared implementation and validation

The physical API and existing natural-unit `iTEBD` interface share
`iLTRGEvolve`; `GetiLTRGM` and `BilayerTraceLTRG` share `iLTRGObserve`.
`GetiLTRGRslt` now handles pure Chi, returns consistent native-grid arrays and
does not create a second interpolated grid. Its natural outputs are distinct
from this physical API; `InterNum` no longer controls that wrapper's grid.
`InitHam_trotter` also assigns right-side Id/Sz2 operators correctly.

Run `tests/test_iltrg_forward.m` from outside the repository with an absolute
artifact directory. It checks free spins, directed DM dimers, transfer gauge
normalization, energy scaling, native/CGS units, natural wrappers and file/memory
restart.

Run `tests/test_iltrg_j1j2.m` with the same external-artifact convention for the
blocked profile. It checks an independently assembled six-spin Hamiltonian
with alternating NN, directed NNN DM and oblique field; grouped free-spin
normalization on all axes; an exact infinite-chain J1z–J2z transfer reference;
energy-scale invariance; the decoupled J1=0 and grouped J2=0 limits; and file/
memory restart with changed-J2 rejection. Its quantum J1–J2 comparisons use
independent eight- and ten-site ED at high temperature, then vary D and tau
separately at common physical temperatures. These are small-D/high-T checks,
not production convergence certification.

Both suites passed in MATLAB R2022b on 2026-11-13. Representative errors from
the J1–J2 run were:

| Independent reference | Maximum discrepancy |
| --- | --- |
| Six-spin physical Hamiltonian vs grouped Hamiltonian | Frobenius norm 9.73e-16 K |
| Classical infinite-chain J1z–J2z transfer matrix | lnZ per spin: 2.34e-15 |
| Quantum J1=1 K, J2=0.4 K vs ten-site ED | lnZ per spin: 5.12e-9 |
| Same quantum case, centered Cm on the same beta grid | 2.46e-6 J/(mol K) |

The quantum comparison used Bz=0.2 T, gz=2, D=24, physical tau=0.005 K^-1
and beta=0.01–0.32 K^-1 (T≥3.125 K). The Cm comparison differentiates the ED
lnZ on the same grid, separating evolution error from derivative discretization.
The suite also compares eight-/ten-site ED, D=16/24 and tau=0.01/0.005 K^-1.
These figures apply to those test points; lower temperatures and stronger
frustration require their own convergence checks. File and memory resumes
reproduced the uninterrupted test trajectory exactly.

Compare D and tau independently at common physical temperatures before interpreting
low-T results. The current time evolution uses the existing first-order Trotter
sequence; this refactor does not upgrade its order.

The grouped path has larger local tensors and a different truncation cut. An
NN-chain choice of D does not establish J1–J2 accuracy at the same D. Check D and
physical tau separately, including representative frustrated parameter points,
before using results in optimization. The J1=0 limit represents two decoupled
J2 chains; their sum must still be divided by the total number of physical spins.

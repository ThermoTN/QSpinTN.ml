# iLTRG forward interface

`[result, state] = iLTRGForward(request, options)` performs an infinite-chain
thermal calculation with a repeated two-site unit cell. The numerical API is
`qspintn.iltrg.forward.v1`, declared in [interface.json](interface.json).
Its currently validated physical profile is spin-1/2 with real nearest-neighbor
bilinear coefficients, including directed DM terms. This does not implement
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

All coefficients and `energy_scale_K` are in K. `tau` is a physical inverse
temperature step in K^-1, so `beta_n=2*n*tau` and `T_K=1/beta_n`. Internally
exchange and field energies are divided by the energy scale and the evolution
step is multiplied by it. `N_max` must cover the target temperature and leave
one extra neighbor for the centered heat-capacity derivative.

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
| `result.summary` | API version, iteration, resume position, actual bond dimensions, elapsed time and numerical identity |
| `result.hamiltonian` | Actual `Hab`, `Hba` and numerical `Para` |
| `result.science` / `result.units` | Request, source hashes and physical conventions |
| second return `state` | Tensor evolution state suitable for `options.state` |

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
restart. The tests are small-D checks, not production convergence certification.
Compare D and tau independently at common physical temperatures before interpreting
low-T results. The current time evolution uses the existing first-order Trotter
sequence; this refactor does not upgrade its order.

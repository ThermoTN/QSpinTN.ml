# tanTRG

`GettanTRGRslt(Para)` builds the Hamiltonian and SETTN initial thermal MPO, then
calls `tanTRGMain(Para,H,rho,Op)`. Existing model, SETTN and temperature-schedule
parameters remain in use. Examples are in `Example_SpinChain_2_tanTRG` and
`Example_TriARX_1_tanTRG`.

```matlab
Para.TanTRGUpdateSites = "cbe"; % 1, 2, "cbe", or "auto"
Para.TanTRGEnvStorage = 'disk';
Para.TanTRGEnvRoot = '/tmp/ygao';
Para.TanTRGKeepEnv = false;
Para.CBERatio = 0.1;
Para.CBEOversampling = 8;
Para.CBEPowerIterations = 1;
Para.CBESvdTol = 1e-12;
Para.CBEExpansionCap = []; % Optional temporary cap; MCrit is the retained cap
Rslt = GettanTRGRslt(Para);
```

`TanTRGUpdateSites="auto"` is the default and preserves the old scheme: start
with two-site evolution; if `Full2=false`, switch to fixed one-site evolution
after the completed sweep exceeds `beta_switch`. Numeric `1` / `2` and `"cbe"`
force the selected scheme from the initial MPO. CBE uses the
[shared RSVD pre-expansion controls](../TensorFunction/CBE.md). Fixed one-site
does not grow the input bonds; CBE and two-site can grow them. The legacy
two-site tensor layout requires `L>=3`.

`TanTRGEnvStorage` defaults to `'memory'`; `'disk'` writes rolling left/right
environments to a unique `tantrg_<job-or-local>_<token>` directory beneath
`TanTRGEnvRoot` (default `/tmp/ygao`). Set a server-appropriate scratch root.
`TanTRGKeepEnv=true` is the default; `false` removes that run's directory on exit,
including errors. `Rslt.EnvDir` reports it in disk mode. Each MAT stores `block`,
`slot`, `side`; writes replace the previous block only after saving succeeds.
The active few environments are loaded outside the Lanczos iterations. SETTN,
the MPO itself and local work arrays remain in memory. These caches are not a
full thermal MPO checkpoint and cannot independently resume an evolution.

`beta` begins at `2*tau0`; every full forward/backward sweep advances it by
`4*tau`. The schedule selected by `GetNowtau` may overshoot `beta_max`, as before.
`rho` represents `exp(-beta*H/2)` with its Frobenius amplitude in `rho.lgnorm`.
`LnZ=2*rho.lgnorm/L`; `En` is energy per site. All entries, including the initial
point, now use the same Hamiltonian scale and per-site normalization.
`EE` and `TrunErr` now align with `beta`, with zero placeholders at initialization;
the former endpoint trimming has been removed. Raw historical results can have
different diagnostic lengths and an unscaled initial energy.

`Rslt.TanTRGInfo.Sweeps` records mode, beta, elapsed sweep seconds, retained bond,
truncation, entropy, added states, expansion seconds and maximum expanded bond.
`TanTRGSweep(rho,H,Para)` runs exactly one full sweep with `Para.tau` and explicit
mode `1`, `2` or `"cbe"`; an optional fourth environment-store argument reuses
prepared environments. Inputs must have the orthogonality center at site 1,
right isometries at sites 2:L, unit Frobenius norm in the tensors, and the global
amplitude in `lgnorm`. This is also the entry point for fixed-input benchmarks.

Optional `GettanTRGCorr` and `savetanTRGRslt` hooks run when present on the path.
Failures in present hooks now propagate, so requested output failures cannot be
mistaken for success. No hook means no automatic full-result file. Result saving
is separate from environment storage. Heat capacity remains a downstream
derivative of energy; compare small systems with ED and check step, retained
rank and initialization convergence before interpreting thermodynamics.

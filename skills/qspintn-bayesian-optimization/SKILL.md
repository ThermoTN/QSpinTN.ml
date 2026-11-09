---
name: qspintn-bayesian-optimization
description: Choose QSpinTN and MATLAB Bayesian optimization settings for exploratory searches, refinement, stalled thermodynamic fits, loss-scale changes and verified continuation.
---

# Choose Bayesian fitting strategies

Use this guidance when selecting or revising a fitting optimizer. Choose a strategy from the evidence; the alternatives below are not a mandatory sequence. A consumer's current tool schema, reviewed plan and numerical-history checks determine what it can execute. Source paths refer to the QSpinTN.ml repository root.

## Checked package interface

[UtilityFunc/OptFunc.m](../../UtilityFunc/OptFunc.m) passes these `QMagenConf.Config.BOPara` fields to MATLAB:

| QSpinTN field | MATLAB option | Meaning in the package |
| --- | --- | --- |
| `AcqFunc_Name{i}` | `AcquisitionFunctionName` | Acquisition function for group i |
| `AcqFunc_ER(i)` | `ExplorationRatio` | Exploration control for a plus acquisition |
| `IsObjDet` | `IsObjectiveDeterministic` | Whether repeated evaluations at fixed inputs are deterministic |
| `Group_MaxEval(i)` | `MaxObjectiveEvaluations` | Initial group limit, then additional evaluations on native resume |
| `Group_Number` | Package loop count | First group calls bayesopt; later groups call resume |
| `Parallel` | `UseParallel` | Serial or package-managed parallel calculation |

[ImportBOPara.m](../../Example/Example_QMagen_SpinChain/ImportBOPara.m) provides examples, not universal defaults: ED has one 50-evaluation group with EI-plus and ratio 0.5; its longer ratio array does not activate unused groups. The iLTRG/XTRG examples have ten groups and decrease the ratio from 0.5 toward 0.1. Ensure group arrays agree with the intended group count and total cost.

The current `OptFunc` external `Restart.Resume` branch contains an interactive `keyboard`; it is not an unattended restart endpoint. It also writes relative `../Tmp/` paths. Use a consumer's controlled adapter or an explicitly prepared native working directory; a skill is not authorization to run the raw example script. No solver implementation is changed by this guide.

## Acquisition and exploration controls

| MATLAB acquisition | Selection emphasis | When to consider it |
| --- | --- | --- |
| `expected-improvement` (EI) | Expected improvement | Refinement while retaining uncertainty-aware search |
| `expected-improvement-plus` | EI plus an overexploitation response | Broader exploration |
| `expected-improvement-per-second` | Time-weighted improvement | Strongly varying forward costs |
| `expected-improvement-per-second-plus` | Time weighting and plus | Cost-aware exploration |
| `lower-confidence-bound` | Posterior mean minus two standard deviations | Alternative mean/uncertainty tradeoff |
| `probability-of-improvement` | Improvement probability with a noise-based margin | Likelihood rather than amount of improvement |

These are the native [acquisition definitions](https://www.mathworks.com/help/stats/bayesian-optimization-algorithm.html). LCB has no public adjustable kappa here; `ExplorationRatio` does not tune LCB or ordinary EI.

**ExplorationRatio controls exploration versus refinement only for plus functions.** Their test is `sigma_F(x) < ExplorationRatio * sigma_noise`: larger values trigger extra exploration more easily; smaller values allow more concentration. It is not a sampling percentage and has no upper bound of 1. Use positive values for portable recipes. Ordinary EI remains uncertainty-aware. See the [plus mechanism](https://www.mathworks.com/help/stats/bayesian-optimization-algorithm.html).

Other controls from the [bayesopt documentation](https://www.mathworks.com/help/stats/bayesopt.html):

- Set `IsObjectiveDeterministic=true` for repeatable evaluations. Experimental noise in fixed data does not itself make f(x) stochastic; assess solver repeatability separately.
- `NumSeedPoints` controls initial random sampling, not a permanent exploration fraction. Balance initial coverage with the forward cost; existing evaluated history is a different input.
- Record RNG state and use a fixed seed for comparisons. Parallel timing and per-second acquisitions compromise reproducibility.
- Evaluation/time limits stop expenditure; exhausting them or a flat short trace does not prove convergence. `MaxTime` may be exceeded by a running evaluation.

Bounds express admissible physics. A reviewed smaller search region can aid refinement but can exclude another basin and is not a confidence interval. MATLAB variable `Transform='log'` changes parameter coordinates, not the objective; real-valued ranges must be strictly positive, so signed exchange/DM intervals cannot use it directly. [Variable documentation](https://www.mathworks.com/help/stats/optimizablevariable.html). QSpinTN's [GetOptVar.m](../../UtilityFunc/GetOptVar.m) currently constructs ordinary bounded variables without exposing Transform; do not invent a consumer option for it.

## Select from the observed situation

| Current evidence, on one unchanged fitting objective | Candidate strategy | Check before choosing |
| --- | --- | --- |
| Little coverage; landscape or competing basins unknown | EI-plus with a larger positive ratio; allocate a finite exploration budget | Coverage, cost and physical bounds; GP uncertainty is not physical parameter uncertainty |
| A useful fit exists; seeking finer parameters | Lower the EI-plus ratio, or use EI; consider log10 total loss | Per-curve residuals and maximum pointwise error, not only a small mean loss |
| Many new distinct evaluations do not improve the historical best | Check numerical/data failures, then compare log10 + stronger refinement against renewed exploration | A cached repeat, failed evaluation, new temperature window or short unlucky batch is not proof of an optimum |
| Good region is narrow but the global landscape is uncertain | Alternate reviewed refinement and exploration batches | Retain full compatible history; report tradeoffs rather than claiming global optimality |
| Error floor or structured residuals remain | Check model, units, temperature sampling and solver accuracy before adding budget | ED should evaluate the selected experimental temperatures directly where supported; TN fidelity needs convergence evidence |

For example, EI-plus ratios `1 -> 0.5 -> 0.1` are a possible decreasing-exploration schedule, not MATLAB-prescribed thresholds. Alternatively, finish with EI, where ratio is ignored. Compare alternatives with the same objective, initial observations and bounded extra budget; report both the historical best and new-batch best. An acquisition change need not change the next candidate.

In autoQMagen, “good” can use the current objective's maximum pointwise relative error strictly below 1%, with no nonzero residual at a zero observation. This is a search target, not proof of physical validity. Its enabled early stop takes precedence; if further precision or landscape sampling is wanted, propose disabling that stop explicitly. A previous high-temperature fit cannot satisfy this test for newly added low-temperature data.

## Log10 plus stronger refinement

For a promising or stalled linear-loss fit, offer `log10(L)` together with reduced-ratio EI-plus or EI as a **reviewable candidate**; it need not wait for a huge dynamic range. First inspect residuals, failed points, reproducibility and bounds. Keep the physical loss definition, data weights and zero handling fixed so the comparison has a clear interpretation.

[LossFunc/loss_func.m](../../LossFunc/loss_func.m) implements QSpinTN `LossConf.Design='log'` after combining curve losses (`'native'` leaves the total unchanged); it does not logarithmically transform the data. Its direct log has no zero floor. autoQMagen uses `loss_mode='log10'` and `log10(max(L,1e-300))`, records raw L and the clamp, and rebuilds the surrogate from rescored observations. For positive L above the floor, ordering and mathematical minimizers are preserved, but the fitted GP and search trajectory change. Negative log scores are normal. Below the floor, scores tie; raw zero loss needs no further loss minimization. Log10 does not repair undefined relative residuals, incorrect units, interpolation errors or an inadequate model.

Examples for autoQMagen's strategy revision, choosing **one** alternative:

```json
{"loss_mode":"log10","optimizer":{"search_mode":"custom","options":{"acquisition_function":"expected-improvement-plus","exploration_ratio":0.1}}}
```

```json
{"loss_mode":"log10","optimizer":{"search_mode":"refine"}}
```

The current autoQMagen `refine` mode maps to ordinary EI, where ratio is ignored; `custom` preserves the specified native options. The ratio 0.1 is illustrative. If the current ratio is already smaller, 0.1 would increase exploration; compare against the active settings. Preserve explicit user choices and the reviewed budget. Do not silently change loss, bounds, stopping policy or solver because these examples exist. Follow an application's existing review gate or applicable prior authorization; do not invent a new approval requirement or re-ask a preference already given. A requested direction such as log10 still needs a concrete plan within that application's workflow.

## Reuse the correct information

Native MATLAB [resume](https://www.mathworks.com/help/stats/bayesianoptimization.resume.html) can change acquisition and ratio. It cannot change `UseParallel`, `NumSeedPoints` or `Initial*` arguments, and its evaluation/time limits specify **additional** work. Native capability does not establish compatibility of a consumer's checkpoints.

| Change | Reusable evidence in autoQMagen's current adapter |
| --- | --- |
| Same scientific objective, parameters and native settings | Verified history; matching native checkpoint resumes with saved RNG |
| Acquisition or ratio changes | All compatible observations seed a new optimizer identity/GP; no repeated physical evaluations needed for those observations |
| Only linear ↔ log10 total loss changes | Rescore saved physical predictions, preserve raw metrics, and rebuild the GP with new scores; old GP scores cannot be resumed unchanged |
| New observations, model, size, method or precision | Transfer suitable parameter candidates and re-evaluate under new conditions; do not import old losses as current observations |

Reusing predictions requires verified, complete matching physical outputs, not merely a table of scalar loss values. The current autoQMagen fitting interface does not expose partial-temperature patching of an old evaluation; when its objective identity changes, use its candidate re-evaluation flow rather than promising only missing temperatures will be computed.

autoQMagen currently exposes acquisition, ratio, deterministic flag and seed-point count as native options, with an application-level seed and batch budget. Its `auto`/`explore`/`refine`/`custom` modes and early stopping are consumer policies. Advanced native settings, variable transformations and parallelism require an implemented adapter, not just a request in a skill.

Explain each recommendation using observed quality, improvement history, objective scale, expected cost and what can be reused. Distinguish saved scientific points from normal optimizer completion; an exit failure is not convergence. Exploration can inform the landscape but does not produce calibrated physical confidence intervals or guarantee escape from every local optimum.

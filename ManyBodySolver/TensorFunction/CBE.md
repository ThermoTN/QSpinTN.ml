# RSVD controlled bond expansion

Select `DMRGUpdateSites="cbe"` or `TanTRGUpdateSites="cbe"`. These are single-site
sweeps with **pre-expansion**, followed by SVD truncation to `MCrit`. The numeric
implementation follows the projected-residual construction in
`tanTRG_Spin/BondExpansionfunc/PreBE_RSVD.m`; it does not depend on QSpace.

For the two-site action, factor the component outside both kept spaces as
`F = (I-U*U') * X * Y * (I-B'*B)`. `U` spans the center's left grouping and `B`
is the right-canonical neighboring tensor. RSVD samples `F' * Omega`, applies
optional power iterations, and takes an SVD of the reduced `F*Q`. The large
two-site matrix `F` is never assembled. Reproject and orthonormalize the selected
right directions, append them to `B`, and zero-pad the center, preserving the
state before the local solve. Reverse the construction for a leftward move.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `CBERatio` | `0.1` | Expanded target `ceil((1+CBERatio)*MCrit)`; `0` disables expansion. |
| `CBEOversampling` | `8` | Extra Gaussian sample columns, bounded by available complements. |
| `CBEPowerIterations` | `1` | RSVD power iterations with QR stabilization. |
| `CBESvdTol` | `1e-12` | Relative singular cutoff for expansion; also the CBE-DMRG truncation cutoff. |
| `CBEExpansionCap` | `[]` | Optional hard cap on the temporary expanded bond; must be at least `MCrit`. |

The added dimension is the target minus the current bond, limited by both
physical complements and nonzero residual singular values. It is not always
`CBERatio*MCrit`, particularly for a small supplied initial state. Full physical
spaces and negligible residuals need no expansion. Complex tensors use complex
Gaussian sketches. `rng` controls reproducibility; the solver does not reset it.

The retained bond after every update is at most `MCrit`; temporary expanded
bonds can exceed it unless `CBEExpansionCap=MCrit` is set. With that cap, bonds
already at `MCrit` cannot expand. CBE-tanTRG uses `TSRGStol` for its post-evolution
SVD, and includes the retained norm and backward bond evolution in `rho.lgnorm`.
Unlike fixed-rank one-site sweeps, CBE has a nonzero truncation error.
The reported discarded weight measures the final SVD, not the residual
directions missed by RSVD or the time-evolution projection error. A tiny
discarded weight alone does not establish thermal time-step accuracy.
`AddedStates` sums added directions across local updates; it is not a maximum
bond dimension.

For speed comparisons, copy the same normalized right-canonical MPS/MPO for
every method, set the same retained cap and local solve/evolution controls,
and time a complete forward/backward sweep. Report input/output/temporary bond
dimensions, truncation and numerical change alongside time. Initialization,
measurements and saving a full result are separate costs; disk sweep timing
includes environment I/O. One sweep is a timing sample, not a convergence or
accuracy-equivalence claim.

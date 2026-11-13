---
name: qspintn-model-building
description: Build and validate QSpinTN spin Hamiltonians and interaction maps for chains, two-dimensional lattices or arbitrary finite graphs, including one-dimensional site ordering, parameter constraints and the distinction between magnetic fields and on-site interactions.
---

# Building a spin model

Source paths below are relative to the QSpinTN.ml repository root. Use this skill when an interaction is absent from an existing model, when changing a finite lattice, or when checking that a displayed Hamiltonian matches the numerical one. A new model name or a skill is not itself an executable solver adapter.

## Define the physical model first

Specify local spin, site numbering, every interaction and bond set, boundary conditions, signs, coupling units and the Zeeman convention. For spin 1/2 use Sα = σα/2; Pauli matrices alone change exchange energies by a factor of four. State whether finite periodic bonds are unique or counted with multiplicity. In particular, the existing XXZtest map counts its L=2 periodic bond twice.

Keep a list of independent parameters, fixed values and equality constraints. Substitute fixed values and resolve ties when displaying the actual fitting model. Remove zero terms, retain fixed nonzero interactions, and show the magnetic field for each dataset. A g component absent from all applied fields is not determined by those data. Do not interpret a fixed parameter as a measured result.

## Choose the model interface

For native QSpinTN models, implement a metadata function `SpinModel/SpinModel_<Name>.m` and an interaction map in `SpinModel/IntrcMap_Folder/`. Review `SpinModel_XXZtest.m`, `IntrcMap_XXZtest.m`, `UtilityFunc/GetModel.m` and `UtilityFunc/PassPara.m` together: metadata and the map must agree on parameter order, spin, geometry, g directions and scale. Declaring `AvlbSolver` does not verify that a solver can execute the model. iLTRG additionally needs its own compatible Trotter construction.

For a consumer that accepts a declarative model, use that consumer's versioned schema and supported operators. autoQMagen retains `aqm.model_spec.v1` for spin-1/2 xx/yy/zz exchange and provides `aqm.model_spec.v2` for all nine real bilinear products `S_i^alpha S_j^beta` (alpha,beta in x,y,z), including DM and symmetric off-diagonal exchange. Both accept chain distances or explicit finite bonds and use ED. It resolves the declaration to a fixed adapter and calls QSpinTN's existing ED path without adding a SpinModel file. Its schema, model hashes, fitting bounds and execution restrictions belong to autoQMagen (`docs/MODEL_SPEC.md` and `autoqmagen/model_spec.py` there); they do not extend QSpinTN's native API. A Heisenberg J1–J2 chain needs three components for each distance, not only SzSz.

## Declarative bilinear models without new package files

When working through autoQMagen, encode a requested bilinear Hamiltonian in `models[].definition` using v2, not in a new `SpinModel_*.m` file. Each term has two operator axes, a real Kelvin parameter, an optional signed `factor`, and a chain or explicit bond set. This represents the general real exchange tensor `sum_ab J_ab S_i^a S_j^b`; parameter names are user-defined, not a model-name whitelist. Bind shared couplings through the same parameter or ties, and supply the requested bounds/fixed values separately.

For directed bonds i→j, expand `D · (S_i × S_j)` as:

| Component | First term | Second term |
| --- | --- | --- |
| Dx | +Dx Sy_i Sz_j | −Dx Sz_i Sy_j |
| Dy | +Dy Sz_i Sx_j | −Dy Sx_i Sz_j |
| Dz | +Dz Sx_i Sy_j | −Dz Sy_i Sx_j |

For example, Dz uses `operators=["x","y"], parameter="Dz", factor=1` and `operators=["y","x"], parameter="Dz", factor=-1` on identical directed bonds. Add xx/yy/zz with the same J for an isotropic Heisenberg interaction. Reversing a bond reverses the DM sign if D is held fixed. Use actual complex Sy matrices for general mixed products; the native real-yy shortcut cannot be applied to a term with a single Sy.

V2 OBC chain links point i→i+r. PBC counts each unordered geometric pair once and keeps the first i→wrapped(i+r) encountered for i=1..L; short periodic samples therefore need particular care. Use explicit directed pairs and separate terms or coefficient factors for a different intended multiplicity. Opposite explicit directions are distinct contributions; never silently deduplicate them as undirected edges. Applied B and g remain in the separate Zeeman path.

A successful declaration still needs forward validation. The consumer's installed tool schema is the executable capability contract; a frozen older knowledge snapshot may describe only v1. Do not reject DM merely because there is no native catalog model, and do not claim unimplemented higher-spin, non-field on-site, multi-spin or TN capabilities.

## Interaction map contract

When generating a map, first construct the physical sites and bonds, then choose a bijection from those sites to `1..L`, and finally emit the automata records. A one-dimensional MPO ordering does not make a two-dimensional or arbitrary graph a nearest-neighbour chain. Essential examples are included below so consumers that snapshot only this SKILL.md still receive the mapping instructions. [Lattice-to-automata mapping](references/mapping.md) provides additional boundary and basis examples for direct repository readers.

### Construct the physical bonds and choose a one-dimensional order

- **1D:** use `p(i)=i`; OBC distance-r bonds are `(i,i+r)` for `i=1..L-r`. For L=4, J1 acts on `(1,2),(2,3),(3,4)` and J2 on `(1,3),(2,4)`. Each Heisenberg bond needs xx, yy and zz components. State the periodic counting rule before adding wrap bonds; unique pairs and periodically translated links can have different multiplicities for short chains.
- **2D:** generate bonds in `(x,y)` coordinates first, then use, for example, `p(x,y)=(x-1)Ly+y`. A square 2×3 OBC lattice has y bonds `(1,2),(2,3),(4,5),(5,6)` and x bonds `(1,4),(2,5),(3,6)`. The latter stay long-range interactions in the MPO. A snake order uses this formula for odd x and `p=(x-1)Ly+(Ly+1-y)` for even x; remap every bond and observable consistently. `IntrcMap_TLXXZ.m` uses y-fast column ordering, with its own triangular offsets and explicit BCX/BCY wrap loops; it is not a square-lattice map.
- **Multiple basis sites:** include `a=1..q`, set `L=q Lx Ly`, and choose `p(x,y,a)=q[(x-1)Ly+(y-1)]+a`. Generate intra-cell and inter-cell interactions with the intended physical/image offsets before applying this numbering. A unit-cell count is not the spin count.
- **Arbitrary lattice or graph:** keep a physical site table and edge list; coordinates are optional. With labels A,B,C,D and edges A–B,A–D,B–C, the order A→1,B→2,C→3,D→4 gives `(1,2),(1,4),(2,3)`. Do not replace the edge `(1,4)` by a path of three short edges. Save both the ordering and inverse ordering so observables and plots can return to physical labels.

For a new finite periodic graph, preserve image information until deciding which duplicate pairs represent multiplicity. Self-edges after wrapping need explicit physical handling, not silent conversion into a two-site term. autoQMagen v1/v2 accepts the resulting pairs via `bonds.kind="explicit"` with `BC="OBC"`: all intended periodic links are already in the list and no additional wrapping is applied. Save the physical coordinate/basis mapping beside the definition.

The current `ManyBodySolver/TensorFunction/AutomataInit.m` consumes a struct array:

```matlab
Intr(n).site = [i, j];             % ascending one-based sites, i < j
Intr(n).operator = {Sx, Sx};       % operators matching the local dimension
Intr(n).J = J_in_natural_units;   % scalar coefficient
```

Convert a physical term `c A(u) B(v)` using the same permutation for sites and operators:

```matlab
sites = [p_u, p_v];
ops = {A, B};
[Intr(n).site, order] = sort(sites);
Intr(n).operator = ops(order);
Intr(n).J = c_natural;
```

For example, `Sz(A) Sx(B)` with `p(A)=5,p(B)=2` becomes `site=[2,5]` and `operator={Sx,Sz}`. Such mixed-axis products are executable through autoQMagen v2; v1 remains restricted to xx/yy/zz. In v2 retain the physical direction in each input pair: `[5,2]` with `["z","x"]` means `Sz(5) Sx(2)`. Its fixed MATLAB adapter sorts both lists together. Operators on the same physical site are local matrix products in their specified order, not a two-site self-edge.

For separated sites, AutomataInit inserts identities between the endpoints. Sort sites together with their associated operators; sorting only the sites changes an interaction when the endpoint operators differ. Each record contributes once; repeated records add their contributions. Longer bonds increase the number of active automata channels across intermediate cuts. Use a whitelist of actual operator matrices when converting external model descriptions; external text must not become an expression passed to `eval`.

The pipeline is `physical graph → fixed site ordering → site/operator/J → AutomataInit → AutomataInit1Site → MPO2H/solver`. For a two-site record, only list its endpoint operators: the constructor propagates identities through intervening sites. Estimate MPO cost from the number of interaction records crossing each cut; alternative orderings can reduce this number while preserving the graph. Reordering does not reduce ED dimension `(2S+1)^L`.

Some older maps return `JmpOut/JmpIn/CS` rather than `site/operator/J`; do not assume they match the current constructor. The checked-in constructor's one-site branch refers to an undefined `JmpIn`. Arbitrary on-site interactions therefore require a separately fixed and tested path. The ordinary Zeeman term uses `AutomataInit1Site.m`, which is a distinct path; do not duplicate it in the interaction map.

The XXZtest map represents y exchange with `Yreal = imag(Sy)` and a negative coefficient: Sy⊗Sy = −Yreal⊗Yreal. This is a paired convention. Copying the real matrices while keeping a positive coefficient reverses the yy coupling. Verify the resulting dense Hamiltonian instead of inferring signs from variable names.

## Energy and field conventions

The legacy direct `GetModel(conf,g,values)` convention selects ES from the absolute
value of `Para_EnScale`, which must be nonzero on that path. K parameters are divided
by ES; parameters marked ES may be dimensionless ratios. Inspect their actual map
use: an XXZ ratio Delta can produce a physical J*Delta term and is not itself K.

A consumer can instead normalize each evaluated candidate using
`E0=max(max(abs(physical interaction coefficients)), max((muB/kB)*norm(g.*B)))`
over all of that candidate's field curves. Expand factors, ties and ratio products
before taking the maximum; keep one E0 across the curves. Only exactly all-zero
interactions and fields use a recorded 1 K reference. autoQMagen's fitting adapter
uses this policy, with no fixed energy_scale_K in the declarative ModelSpec. It
must preserve the physical map and rescale every energy term, not merely replace
conf.ES. Numerical convergence and supported-model checks remain independent.
A change of reference alone preserves physics only when physical temperatures,
fields and solver accuracy (including physical imaginary-time steps) are retained.

`PassPara` supplies `Para.Model`, `Para.Geo`, `Para.d`, `Para.L`, `Para.IntrcMap_Name` and the field/g information. Physical input B is in tesla. The existing field path yields HZ/kB = −(μB/kB)Σi,α gα Bα Siα. Read [units](../qspintn-units/SKILL.md) before converting thermal results, and [ED](../qspintn-ed/SKILL.md) for the current finite-field M/H and output-grid conventions.

Treat single-ion terms such as `D(Sz)^2` and a model coefficient `−Delta Sx` separately from a measured applied magnetic field. They use model energy units and local operators, not automatic tesla/g conversion. An explicitly intended transverse magnetic field can instead use the Bx path; never include both descriptions for the same physical contribution. For spin 1/2, `(Sz)^2 = I/4`, so uniform D contributes only a constant and cannot be inferred from Cm or finite-field χ. Arbitrary non-field on-site terms are currently outside autoQMagen's v1/v2 adapters and encounter the native one-site constructor issue described above; a skill does not make that numerical path available.

## Evidence before fitting

Use the smallest informative system and independently assemble the dense matrix from tensor products. Check Hermiticity, spectra, boundary-bond multiplicity and a nonzero-field case. Compare a limit that reproduces an existing model, then a point where the new interaction is nonzero. Test zero couplings if the chosen scale convention permits them. Compare physical Cm and χ, including per-spin molar normalization; interpolation errors need tolerances separate from spectral roundoff.

Only then run a small synthetic parameter-recovery test through the consumer's normal optimizer and report path. A successful software recovery does not establish a material model, parameter identifiability or finite-size convergence. Save the definition, actual bond list, solver settings, source revision/hashes and validation evidence with the run. Changes to any of these require a new numerical identity rather than reuse of old cached results.

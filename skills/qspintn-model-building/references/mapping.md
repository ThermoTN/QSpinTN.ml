# From physical sites to QSpinTN automata records

Read this reference when generating or reviewing an interaction map. All source paths are relative to the QSpinTN.ml root. The examples define finite Hamiltonians; they do not infer interactions or lattice symmetry from a material name.

## 1. Keep geometry, ordering and operators separate

Retain three explicit objects:

1. **Physical sites:** stable labels, coordinates if known, basis/site type, and local spin/basis. A graph does not require Euclidean coordinates.
2. **Physical interactions:** endpoint labels, interaction class, actual endpoint operators, coupling, sign, and multiplicity. Include the intended periodic image displacement where relevant.
3. **Ordering:** a bijection `p(physical_site) -> 1..L` and its inverse. Use the same mapping for observables, field directions, output labels and validation.

Build interactions in physical coordinates or labels before applying the ordering. A geometric near-neighbour can become a long-range one-dimensional MPO bond; preserve it. Distances in the ordered list do not define physical neighbours in a two-dimensional lattice or irregular graph.

For a pair `c A(u) B(v)`, translate to `i=p(u)`, `j=p(v)` and reorder the associated operators when needed:

```matlab
sites = [i,j];
ops = {A,B};
[ordered_sites, permutation] = sort(sites);
Intr(n).site = ordered_sites;
Intr(n).operator = ops(permutation);
Intr(n).J = c_natural;
```

For example, `c Sz(A) Sx(B)` with `p(A)=5`, `p(B)=2` becomes `site=[2,5]`, `operator={Sx,Sz}`. It does not become `{Sz,Sx}` on `[2,5]`. Distinct-site spin operators commute, so this paired permutation preserves the physical term. Directed antisymmetric couplings need their full signed operator expansion before canonicalization; treating them as undirected Heisenberg bonds changes the model. autoQMagen v1 rejects such mixed operators; v2 accepts all nine xyz pair products and retains input bond orientation until the adapter sorts sites and operators together.

Two operators on the same physical site must be multiplied in their intended order into a local matrix; they are not a two-site self-edge. Whether that local matrix can be executed depends on the on-site path described below.

## 2. One-dimensional chain

For L sites use `p(i)=i`. With open boundaries, a distance-r pair has `i=1..L-r`, `j=i+r`. A Heisenberg coupling on a bond is three records, `J Sx_i Sx_j`, `J Sy_i Sy_j`, and `J Sz_i Sz_j`. An XXZ coupling gives the first two coefficient Jxy and the third Jz.

For L=4 OBC, a J1–J2 chain has:

| Class | Ordered pairs | Components |
| --- | --- | --- |
| J1 | (1,2), (2,3), (3,4) | xx, yy, zz |
| J2 | (1,3), (2,4) | xx, yy, zz |

The J2 interactions must not disappear merely because their sites are not adjacent in the MPO.

For periodic chains, state the counting convention before generating wrap bonds. A graph of unique unordered pairs and a sum over periodically translated directed lattice links can differ for very small L. Native `IntrcMap_XXZtest.m` appends `(1,L)` to the OBC list: L=2 PBC therefore counts `(1,2)` twice. autoQMagen's custom chain declaration uses unique unordered pairs; its L=2 PBC model is deliberately not the native double-counted convention. To represent a physical multiplicity in a unique-pair schema, use an explicit coefficient factor and record why.

Native XXZtest is a useful operator/sign reference, but its final loop uses `length(int_cell)` rather than the number of rows. In the checked-in snapshot L=2 OBC has three rows and five columns, so that implementation indexes beyond its rows. Do not copy this loop bug into a new map; iterate the actual record count (`size(int_cell,1)` or the filled count).

## 3. Two-dimensional lattice

For one site per cell, a convenient column-by-column ordering is

\[
p(x,y)=(x-1)L_y+y,
\qquad 1\le x\le L_x,\quad1\le y\le L_y.
\]

This is the ordering used by `IntrcMap_TLXXZ.m`: y varies fastest. For Lx=2, Ly=3, the site-to-order table is:

| Physical site | Column ordering | Snake ordering |
| --- | ---: | ---: |
| (1,1) | 1 | 1 |
| (1,2) | 2 | 2 |
| (1,3) | 3 | 3 |
| (2,1) | 4 | 6 |
| (2,2) | 5 | 5 |
| (2,3) | 6 | 4 |

The alternate snake ordering is `p=(x−1)Ly+y` for odd x and `p=(x−1)Ly+(Ly+1−y)` for even x. Choosing it requires remapping every physical interaction and observable; it is not compatible with reusing an unmodified map that assumes the column ordering.

For a square 2×3 OBC nearest-neighbour graph, the column ordering yields:

- y edges: `(1,2), (2,3), (4,5), (5,6)`;
- x edges: `(1,4), (2,5), (3,6)`.

The same graph in snake order has y edges `(1,2), (2,3), (5,6), (4,5)` and x edges `(1,6), (2,5), (3,4)`. Both contain seven physical bonds and have identical spectra under a site permutation, although their MPO construction costs can differ.

For the triangular map, consult the actual coordinate offsets rather than copying the square graph. Its OBC body uses first-class offsets `(0,1), (1,0), (1,1)` and second-class offsets `(1,−1), (1,2), (2,1)` when the destination remains inside the sample. On the 2×3 column ordering, the first-class diagonal edges add `(1,5), (2,6)`; the second class has `(2,4), (3,5), (1,6)`. These offset classes describe the checked-in map, not a universal choice of triangular primitive vectors.

Periodic boundaries require the intended cell identifications as well as the offsets. Generate wrap bonds from the physical lattice, retain their image/multiplicity information until the counting rule is applied, and then map to ordered sites. Small circumferences can merge two image bonds onto the same site pair or create self-edges. Decide physically whether to sum, retain multiplicity, or reject those cases; unconditional deduplication changes some finite periodic Hamiltonians. Native TLXXZ has explicit BCX/BCY supplement loops; reproduce and check those loops when claiming equivalence, rather than assuming generic modulo wrapping produces the same finite graph.

For q sites per unit cell, include a basis label `a=1..q`, set `L=q Lx Ly`, and for example use `p(x,y,a)=q[(x−1)Ly+(y−1)]+a`. Build intra-cell and inter-cell edges using the basis offsets, then apply boundary identifications. A unit-cell count is not a spin count. Inspect each native map's implemented boundary fields: for example, a declared BCX field alone does not establish that its map adds x-periodic edges.

For autoQMagen v1/v2, freeze the resulting finite pairs as `bonds.kind="explicit"` and use the finite `L` with `BC="OBC"`; periodic edges, if any, must already appear in that explicit list. Here OBC means no additional automatic wrapping. Save the physical coordinate/basis table beside the declaration, because the numerical schema contains the ordered graph, not crystallographic coordinates.

## 4. Arbitrary finite graph and automata cost

For labels `{A,B,C,D}` with physical edges `{A–B,A–D,B–C}`, choose for example `p(A)=1,p(B)=2,p(C)=3,p(D)=4`. The map contains `(1,2),(1,4),(2,3)` for each required operator component. The `(1,4)` interaction propagates identities through sites 2 and 3; it is not approximated by edges 1–2–3–4. Site coordinates, irregular bond lengths and missing links may be retained as metadata without changing this representation.

The active implementation performs:

```text
physical interaction graph
    -> fixed site ordering and associated operators
    -> Intr(n).site / operator / J
    -> AutomataInit (identity propagation between interaction sites)
    -> AutomataInit1Site (external field)
    -> ED MPO2H or a compatible tensor-network solver
```

`AutomataInit.m` increases channel counts on the first/last interaction sites and every cut between them. In this direct constructor, each pair record crossing a cut contributes a channel in addition to the base channels; three separately stored spin components can contribute three. Long-range edges and wide two-dimensional strips can therefore make a large MPO even when each site has few physical neighbours. Compare possible orderings by the number of interactions crossing each cut; do not remove physical bonds to reduce that count. Reordering does not reduce the ED Hilbert dimension `(2S+1)^L`. Whether a downstream tensor-network implementation compresses efficiently or converges requires separate validation.

The modern automata interface accepts records with more than two distinct ordered sites, but a consumer may expose a narrower schema. autoQMagen v1 accepts spin-half diagonal pair exchange; v2 also accepts general real bilinear exchange and DM on directed bonds. Neither version enables multi-site interactions or general local dimensions.

## 5. External magnetic field and other on-site terms

An applied magnetic field enters the established package path as physical `Field.B` in tesla. `PassPara.m` calls `GetField.m`, which constructs the natural field magnitude using `GetUnitCon.m`, and passes the g factors/directions. `AutomataInit1Site.m` adds the Zeeman operator and constructs the magnetization operator:

\[
H_Z/k_B=-\frac{\mu_B}{k_B}\sum_{i,\alpha} g_\alpha B_\alpha S_i^\alpha.
\]

For each dataset preserve Bx, By, Bz and the measured direction. Do not add the same Zeeman term again in Intr. The displayed formula should substitute the actual field and fixed g values. A g component that never enters any applied field does not become an independently constrained quantity merely because it appears in generic metadata.

Other local physics is a separate model contribution:

| Term | Parameter convention | Implementation distinction |
| --- | --- | --- |
| `D (Sz)^2` | D in K, or explicitly declared natural energy units | Build the local product matrix; it is not a magnetic Bz field |
| `−Delta Sx` | Delta in model energy units | Independent model transverse term unless explicitly identified with an applied field |
| Physical transverse magnetic field | Bx in T and gx dimensionless | Use the external-field path; an equivalent energy coefficient is `(μB/kB) gx Bx` |

For S=1/2, `(Sz)^2=I/4`; uniform D adds `LD/4` to every energy. This affects the absolute energy origin but neither Cm nor finite-field M/H for this spin-only model, so those data cannot determine D. For S≥1 the term can be nontrivial, requiring the appropriate local basis and supported solver path.

The checked-in `AutomataInit.m` one-site branch adds to `H.A{JmpIn}(1,end,:,:)` although `JmpIn` is undefined there. It is therefore not a usable generic non-field on-site implementation. `AutomataInit1Site.m` is a different code path and supports the existing Zeeman handling. To enable arbitrary local terms, first implement or repair a dedicated numerical path within the authorized project, then test a minimal system, its matrix and its response to field. autoQMagen v1/v2 deliberately rejects such declarations; the current skill describes the distinction and required work without claiming it is already supported.

## 6. Check the generated map before fitting

Save the site/inverse-site table and the explicit interaction list with their units and counting rule. Check that every intended physical edge maps exactly once with the intended multiplicity, every operator stays attached to its physical site, and no out-of-range or accidental self-site record is emitted.

On a small sample, assemble an independent dense Hamiltonian from tensor products. Check Hermiticity and eigenvalues, including nonzero fields and the new coupling. Reorder the sites and verify spectral invariance. Compare physical Cm and χ against a known-model limit with the same boundary counting and g/energy conventions. The real-MPO `yy` convention must be checked as a pair: `Yreal=imag(Sy)` requires a negative coefficient because `Sy⊗Sy=−Yreal⊗Yreal`.

Source entry points reviewed for this guidance: `SpinModel/IntrcMap_Folder/IntrcMap_XXZtest.m`, `SpinModel/IntrcMap_Folder/IntrcMap_TLXXZ.m`, `UtilityFunc/GetModel.m`, `UtilityFunc/PassPara.m`, `UtilityFunc/GetField.m`, `UtilityFunc/GetUnitCon.m`, `ManyBodySolver/TensorFunction/AutomataInit.m`, `ManyBodySolver/TensorFunction/AutomataInit1Site.m`, and `ManyBodySolver/ED/GetEDRslt.m`. Recheck behavior after changing the dependency revision.

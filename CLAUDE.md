# spaceTime4foam: instructions for Claude Code

## 1) What this repository is

`spaceTime4foam` is a small OpenFOAM-based toolbox for **space-time finite
volume methods**. Time is treated as an extra coordinate, so an unsteady
problem becomes a *steady* problem on a mesh of one higher dimension.

The name deliberately mirrors `solids4foam`, and so does the repository
style (see section 9).

**Step one (the only task for now):** a 1D+t benchmark comparing

- cell-centred finite volume (CCFV) in space-time, and
- vertex-centred (node-centred, edge-based) finite volume (VCFV) in
  space-time,

against an analytical solution, measuring observed orders of accuracy under
mesh refinement and comparing them with the theoretical orders. Step one
ends with a from-scratch explanatory guide (section 8).

Do not start later steps (diffusion, periodic-in-time, 2D+t, adaptation)
unless the developer asks.

## 2) About the developer

CFD PhD researcher, experienced with OpenFOAM (cell-centred FV, PIMPLE,
dynamic meshes, custom boundary conditions, 3D-0D coupling). Explain design
choices tersely in conversation. Ask before major architectural decisions.

## 3) First actions in a new environment

1. Detect the OpenFOAM fork and version (`echo $WM_PROJECT
   $WM_PROJECT_VERSION`) and report it before writing code.
2. Check that Gmsh, gnuplot and `markdownlint` are available; report what
   is missing.
3. Read the papers in `refs/` (section 10) and the solids4foam guidelines
   file if present in `docs/` (section 9).
4. Propose a short plan for step one and wait for approval.

Check real headers in `$FOAM_SRC` rather than guessing API signatures.

## 4) Naming convention

- Repository: `spaceTime4foam`
- Solver application: `spaceTime4Foam` (compare `solids4Foam`)
- Model library: `src/spaceTime4FoamModels`
- Runtime-selectable base class: `spaceTimeModel`, with derived types
  `cellCentred` and `vertexCentred`, selected in
  `constant/spaceTimeProperties`, so the same case runs with either method
  by changing one keyword
- Runtime-selectable `analyticalSolution` base class (first derived type:
  `travellingSine`; later others), used for boundary data, initial data and
  error norms by both methods
- Geometry class for VCFV: `medianDualMesh` (do not use the name
  `edgeMesh`, which already exists in OpenFOAM)
- Utilities: `perturbSpaceTimeMesh` (random interior node perturbation) and
  `spaceTimeErrors` (error norms against the analytical solution)
- Boundary conditions for CCFV: `spaceTimeAnalyticalFixedValue`, which
  evaluates the selected `analyticalSolution` on a patch, and
  `spaceTimeLinearExtrapolation` (outflow variant, section 7)

Suggested layout, following solids4foam:

```text
spaceTime4foam/
  Allwmake  Allwclean  Alltest
  applications/solvers/spaceTime4Foam/
  applications/utilities/perturbSpaceTimeMesh/
  applications/utilities/spaceTimeErrors/
  applications/test/Test-<name>/
  src/spaceTime4FoamModels/
    spaceTimeModels/spaceTimeModel/
    spaceTimeModels/cellCentred/
    spaceTimeModels/vertexCentred/
    medianDualMesh/
    analyticalSolutions/
    boundaryConditions/spaceTimeAnalyticalFixedValue/
  tests/<testName>/
  tutorials/advection1D/travellingSine/
  docs/
  docs/audits/
  refs/
```

Tests: test applications live in `applications/test/Test-<name>/`
(OpenFOAM.com convention, executable `Test-<name>`). Test cases live in
`tests/<testName>/`, each with an `Allrun` that exits non-zero on failure.
`Alltest` builds nothing itself; it runs every test case and reports a
pass/fail summary. The PDFs in `refs/` are not tracked by git;
`refs/README.md` gives the full citations.

## 5) Space-time conventions

- Mesh x is physical x; mesh y is time t; mesh z is one cell thick with
  `empty` front and back patches (a standard OpenFOAM 2D case).
- Work nondimensionally, so treating t as a length passes dimension checks.
- Patch names: `tStart` (t = 0), `tEnd` (t = T), `xLeft`, `xRight`,
  `frontAndBack`.
- The solver is steady: `runTime` counts iterations, not physical time.

## 6) Step one: the 1D+t benchmark

### Model problem

Linear advection, `du/dt + a du/dx = 0`, on x in [0, 1], t in [0, 1],
with a = 1.

Analytical solution (method of characteristics):
`u(x, t) = sin(2 pi (x - a t))`. Initial data at t = 0 and inflow data at
x = 0 are both taken from this expression, so the data are consistent.

In space-time this is the steady problem `div_st(A u) = 0` with space-time
velocity `A = (a, 1)`.

Boundary treatment (primary runs): data only on inflow boundaries.

- `tStart` and `xLeft`: analytical values (inflow).
- `tEnd` and `xRight`: no data (outflow).

### Schemes to compare

- **CC-1:** CCFV, `upwind`. Theoretical order 1.
- **CC-2:** CCFV, `linearUpwind grad(u)` with `grad(u)` from
  `Gauss linear`. Theoretical order 2. This is the primary CC-2.
- **CC-2-LS:** as CC-2, but with `grad(u)` from `leastSquares`. Extra
  variant, run on both mesh families.
- **CC-2-EX** and **CC-2-LS-EX:** as CC-2 and CC-2-LS, but with the
  `spaceTimeLinearExtrapolation` outflow condition (section 7) on `tEnd`
  and `xRight` instead of `zeroGradient`. Added in S2b, both with the
  corner fallback (section 7). CC-2-EX (Gauss) removes the boundary-strip
  signature. CC-2-LS-EX is kept but is **not** a second-order outflow
  variant: it retains the first-order outflow strip (section 7), and
  must be reported as such.
- **VC-1:** VCFV, upwind flux without reconstruction. Theoretical order 1.
- **VC-2:** VCFV, upwind flux with linear reconstruction
  `uL = uj + 0.5 grad(u)_j . (xk - xj)`. Theoretical order 2.

### Meshes

- Structured right triangles from Gmsh (transfinite, one diagonal
  direction), N x N squares, N = 8, 16, 32, 64, 128 (256 if affordable),
  extruded one layer and imported with `gmshToFoam`.
- Benchmark families use **`Left` diagonals** (from upper-left to
  lower-right, direction (1, -1)). Reason: with `Right` diagonals
  (direction (1, 1)) and a = 1, every diagonal is parallel to A, so
  `A . n = 0` on every diagonal face and face values are copied exactly
  along characteristics (confirmed in the S2 audit). That set-up
  unfairly favours CCFV.
- The `Right` family is kept only as a small aligned case in S5, labelled
  as a characteristic-alignment illustration, never as a main result.
- The same meshes perturbed by `perturbSpaceTimeMesh`: interior nodes moved
  randomly by at most 0.2 h in x and t, boundary nodes not moved (as in
  Tufillaro et al.), front/back node pairs moved identically, fixed random
  seed.
- Both methods must run on exactly the same triangle meshes.

### Error measures (`spaceTimeErrors`)

- L1, L2 and Linf errors over the whole space-time domain, evaluated at each
  method's own unknown locations and weighted by its own control volume
  (cell area for CCFV, dual area V_j for VCFV).
- The same norms reported separately for interior and boundary unknowns,
  to show whether the boundary closure limits accuracy. CCFV: boundary
  cells are cells with a face on a non-empty patch. VCFV: boundary nodes
  are nodes on any boundary edge, including the t = T boundary.
- L1 and L2 are normalised by the area of their own subset
  (`L1 = sum |e| A / sum A`, `L2 = sqrt(sum e^2 A / sum A)`); Linf is a
  plain max. Mesh size `h = sqrt(A_domain / (nTriangles / 2))`, which is
  1/N on both mesh families.
- Order checks are report-only until regression values are frozen after
  S5. Never tune numerics to hit a target order; investigate and report
  unexpected orders instead.
- Final-time error at t = T. For VCFV the nodes lie on t = T. For CCFV,
  document how the face value is obtained (for example, linear
  extrapolation with the cell gradient) and why this matters.
- Observed order between successive meshes:
  `p = log(e_coarse / e_fine) / log(h_coarse / h_fine)`.
- Report errors against h and against the number of unknowns
  (about 2 N^2 cells for CCFV versus (N + 1)^2 nodes for VCFV).
- Also record iterations to convergence and wall time.

Write all results as plain-text `.dat` files and produce convergence plots
(PNG) with gnuplot. Ask before using Python for anything.

### Verification tests (must pass before benchmarking)

1. **CCFV sanity check.** On a quad (blockMesh) space-time mesh, CC-1 must
   match stock `scalarTransportFoam` with `Euler` time stepping on the
   matching 1D mesh to solver tolerance (space-time CCFV with upwinding in
   time is backward Euler plus upwind in space).
2. **Geometry unit test.** Single right triangle with nodes (0, 0), (1, 0),
   (0, 1): every node gets V = 1/6, and the directed-area vector of the
   edge from (0, 0) to (1, 0) is (1/3, 1/6).
3. **Geometry identities** on structured and perturbed meshes, to
   round-off:
   - sum of V_j equals the domain area;
   - closure at every node:
     `sum_k n_jk + 0.5 * sum_(boundary edges at j) n_B = 0`;
   - interior nodes: `V_j = 0.25 * sum_k (p_k - p_j) . n_jk`.
4. **Freestream preservation:** constant u gives zero VCFV residual.
5. **Linear exactness:** linear u gives exact LSQ gradients and an exact
   interior flux divergence.
6. **Published reference:** reproduce the 2D manufactured-solution
   results of Tufillaro et al. (their Table 2, strong Dirichlet on all
   boundaries, `u1 = x^2 + xy + y^2` etc.) with VC-2.

### Done criteria for step one

- `./Allwmake` builds cleanly; `./Alltest` runs all verification tests.
- `tutorials/advection1D/travellingSine/Allrun` runs the full refinement
  sweep for all schemes in section 6 on both `Left` mesh families (plus
  the small `Right` aligned case) and produces the tables and plots.
- The guide in section 8 is written, lint-clean and uses only numbers
  produced by these scripts.

## 7) Numerics reference

### CCFV (`cellCentred`)

Solve `fvm::div(phiST, u) == 0` with `phiST = fvc::flux(Ust)`,
`Ust = (a, 1, 0)`. Use a nonsymmetric linear solver (`PBiCGStab` with
`DILU`) and iterate until converged (deferred correction for
`linearUpwind`). Inflow patches use `spaceTimeAnalyticalFixedValue`;
outflow patches use `zeroGradient` (primary, as specified).

Known effect (S2 audit): OpenFOAM's `linearUpwind` adds no correction on
non-coupled boundary faces, so with `zeroGradient` the outflow face value
is `u_f = u_P`, which is O(h) wrong in the strip of outflow boundary
cells. This gives the boundary-strip signature L1 ~ h^2, L2 ~ h^1.5,
Linf ~ h. OpenFOAM's `leastSquares` gradient is also inconsistent at
non-coupled boundary faces (it uses only the patch-normal part of the
face offset), which makes CC-2-LS first order in inflow boundary cells.
Both effects are reported, not hidden.

Outflow variant (before S5): `spaceTimeLinearExtrapolation`, a boundary
condition setting `u_f = u_P + grad(u)_P . (C_f - C_P)` with the full
offset vector (not only its normal part), re-evaluated every outer
(deferred-correction) iteration from the current `grad(u)`. Its update
mechanism must be documented and its effect on convergence (iterations,
final residual) reported against `zeroGradient`.

Corner fallback (S2b decision): in any cell with two or more
`spaceTimeLinearExtrapolation` faces (on `Left` meshes only the top-right
outflow corner cell), those faces use the zeroGradient value `u_f = u_P`
instead. Reason (S2b audit): in that cell the extrapolation feeds the
(1, -1) gradient component back into itself, a direction the a = 1 cell
balance cannot see. With `Gauss linear` the gain is exactly 1 (neutral
mode, non-unique solution); with `leastSquares` it is 1.5 (divergence
hidden by the residual). The fallback applies to both CC-2-EX and
CC-2-LS-EX. Each run logs the number of faces using the fallback. A
`FatalError` guard stops any configuration in which the unstable corner
coupling could occur without the fallback. Cost: one O(h) cell when
a != 1, which would then dominate Linf; for a = 1 the corner cell value
is unaffected. Upgrade path for a != 1 (for example a multi-rate
benchmark): take the corner cell's gradient from its upwind neighbour.

With the fallback, CC-2-LS-EX converges stably (contraction about 0.91,
independent of N, from a local mode near the outflow corner), but it is
still first order in the outflow strip (S2b audit). Cause: v2412
`leastSquares` uses the normal-only boundary offset `fvPatch::delta()`,
while the extrapolation uses the full offset, so the combination is not
linearly exact and the tEnd/xRight-strip gradient error is O(1). This is
an outflow effect, distinct from the `leastSquares` inflow inconsistency
above. Do not tune it; report it.

Open (S2b audit, minor): the `leastSquares` guard detects schemes by
name; an allow-list of Gauss-type schemes (or a type check) would be
more robust.

The `cellCentred` source term: until a tested source implementation
exists, `cellCentred` must stop with a fatal error if the selected
`analyticalSolution` has a non-zero source.

### VCFV (`vertexCentred`)

Topology: take the triangulated `front` patch. Use `localPoints()`,
`localFaces()`, `edges()` (internal edges first, see `nInternalEdges()`)
and `faceEdges()`. Nodes are patch points; edges are patch edges; store one
vector per edge oriented from the smaller to the larger node label.

Geometry (d = 2):

- Triangle area |T|; dual area `V_j = (1/3) sum_(T contains j) |T|`.
- `n_j^T`: normal of the triangle edge opposite node j, length equal to
  that edge's length, oriented away from j.
- Interior contribution: `n_jk = (1/3) sum_(T contains jk) n_j^T`.
- Boundary correction: for each boundary edge, add `n_B / 6`, where `n_B`
  is its outward normal with length equal to the edge length. This was
  checked by hand on a single right triangle; the closure test must still
  confirm it on full meshes.
- Orientation: determine every normal's sign explicitly with a dot-product
  test; do not rely on Gmsh face winding.

Gradients: unweighted linear least squares over edge neighbours. Store each
node's 2x2 matrix as a `symmTensor` with zz = 1 so `inv()` can be used.

Residual (single edge loop), with `A = (a, 1)` and `nHat = n_jk / |n_jk|`:

```text
Phi_jk = 0.5 (nHat . A)(uL + uR) - 0.5 |nHat . A| (uR - uL)
Res_j += Phi_jk |n_jk|        Res_k -= Phi_jk |n_jk|
```

Boundary closure (weak): each boundary edge contributes to each of its two
nodes with area `|n_B| / 2`, using an upwind flux between the nodal value
and a boundary state (analytical value on inflow, the nodal value on
outflow). Step one uses this simple nodal closure only. The linearly
exact boundary flux quadrature referenced by Tufillaro et al. (Nishikawa,
J. Comput. Phys. 281, 2015) is deferred until that paper is added to
`refs/`; its weights must then be verified from the paper.

Solver: two-stage Runge-Kutta pseudo-time as in Tufillaro et al.
(CFL = 0.5, `dTau_j = CFL V_j / sum_k 0.5 |n_jk . A|`). If convergence is
too slow on fine meshes, propose defect correction with a first-order
Jacobian in an `lduMatrix` built on an `lduPrimitiveMesh`.

Output: write u as a `pointScalarField` (front and back points share (x, t);
match them by coordinates) so ParaView shows the result directly.

## 8) The explanatory guide (deliverable)

Write `docs/spaceTime4foamGuide.md`. Audience: a CFD researcher who knows
OpenFOAM but has never seen space-time or vertex-centred methods. Explain
from scratch, define every symbol, and prefer small worked examples.

Required sections:

1. Overview: what space-time methods are and why this benchmark exists.
2. The model problem and its analytical solution, derived by the method of
   characteristics, including why the inflow data must be consistent with
   the initial data.
3. Space-time reformulation: the space-time flux, inflow and outflow
   boundaries from the sign of `A . n`, and why the initial condition is
   just an inflow boundary condition.
4. CCFV in space-time: integration over a space-time cell, face fluxes,
   upwind and `linearUpwind`, the exact equivalence with backward Euler on
   quad meshes, and how the OpenFOAM implementation maps onto this.
5. VCFV: median-dual control volumes, the dual-area and directed-area
   formulas, the worked right-triangle example, the boundary correction,
   the closure identity, least-squares gradients, the edge-loop residual,
   boundary closure and the pseudo-time solver.
6. From 1D+t to 2D+t:
   - CCFV: a 3D OpenFOAM mesh with z = t and space-time velocity
     (u_x, u_y, 1).
   - VCFV: tetrahedra, `V_j = (1/4) sum |T|`,
     `n_jk = (1/6) sum n_j^T`, an expected boundary coefficient of 1/12
     that must be verified by the closure test, and the time-aligned-edge
     issue for diffusion (Padway and Nishikawa).
   - A short note on 3D+t and why node-centred schemes are favoured there
     (roughly 20 pentatopes per node).
7. Code walkthrough from scratch: every class and file, how data flows from
   mesh to geometry to residual to solution to errors, and exact build and
   run commands.
8. Verification results: each test in section 6 and its outcome.
9. Benchmark results: tables and plots of errors, observed versus
   theoretical orders, accuracy per unknown, iterations and wall time.
10. Findings and nuances: fairness of the comparison (unknown counts,
    evaluation locations, final-time values for CCFV), boundary treatment
    effects, perturbed versus structured meshes, anything surprising.
    Must include the characteristic-alignment case (`Right` diagonals
    with a = 1), the boundary-strip error signature (L1 ~ h^2,
    L2 ~ h^1.5, Linf ~ h) compared with the observed orders, the
    `leastSquares` boundary inconsistency, and what the extrapolated
    t = T value does and does not measure. Also the outflow-corner
    mode of `spaceTimeLinearExtrapolation` (neutral with Gauss, gain
    1.5 with leastSquares), the corner fallback, its O(h) cost when
    a != 1, and the upwind-neighbour-gradient upgrade path; and why
    CC-2-LS-EX keeps a first-order outflow strip (normal-only
    `leastSquares` boundary offset against full-vector extrapolation)
    while CC-2-EX does not.
11. When each approach is likely to be useful, and recommended next steps.
12. References.

Rules for the guide:

- Every number must come from a script in this repository; state the
  command that regenerates it. Never invent or estimate results.
- If a result contradicts theory, report it honestly and investigate.
- Mathematics may use plain-text formulas or LaTeX in `$$` blocks.
- The file must pass `markdownlint`.
- Explain the mesh z thickness of 1: cell volumes equal triangle areas,
  so the CCFV (cell area) and VCFV (dual area V_j) error weights are
  directly comparable.

## 9) Coding style: follow solids4foam

The developer uses the solids4foam contribution guidelines as the style
reference. If the full guidelines file is in `docs/`, follow it. Key points:

- Follow OpenFOAM/solids4foam C++ style: indentation, braces, comment
  style, naming. Keep code readable and explicit; avoid clever or
  condensed expressions and abstraction-heavy designs.
- Comment only non-obvious behaviour.
- Every C++ file carries a standard OpenFOAM-style licence/header block,
  consistent across the repository; keep include ordering consistent.
- Runtime selection: `TypeName("...")` in headers,
  `addToRunTimeSelectionTable(...)` in sources, source files added to
  `Make/files*`, and dictionary `type` strings matching `TypeName` exactly.
- Keep behaviour dictionary-driven and runtime-configurable.
- Build lists follow the solids4foam pattern: `Make/files.openfoam` and
  `Make/files.foamextend`. The detected fork is the primary target; guard
  fork-specific API differences the way solids4foam does, and never claim
  an untested fork works.
- Shell scripts (`Allwmake`, `Allrun`, `Alltest`, `Allclean`) follow
  standard OpenFOAM/solids4foam script style.
- Markdown is concise, practical and passes `markdownlint`.
- Make the smallest correct change for each task; do not reformat or
  rename unrelated code.
- Deliver work as small, focused, reviewable commits, each exportable as a
  git patch. Do not silently change numerics; flag any deviation from the
  papers.

## 10) References (PDFs in `refs/`)

- Tufillaro, Williams and Nishikawa, "Edge-Based Discretizations on
  Triangulations in Any Dimension, With Special Attention to
  Four-Dimensional Space", Int. J. Numer. Meth. Fluids (2026),
  doi:10.1002/fld.70084. Median-dual geometry, directed-area vectors,
  scalar advection verification in 2D, 3D and 4D.
- Padway and Nishikawa, "An Adaptive Space-Time Edge-Based Solver for
  Two-Dimensional Unsteady Viscous Flows", AIAA SciTech 2022. Space-time
  flux splitting, time-aligned-edge fix, pseudo-time solver.

## 11) Agent workflow

The main session is the **orchestrator**. It plans, delegates, tracks
progress and reports to the developer. It does not write large amounts of
code itself. Three project subagents live in `.claude/agents/`:

- `foam-builder`: implements infrastructure and C++/OpenFOAM code.
- `foam-auditor`: independent, read-only review of correctness, numerics,
  verification and style. It reports; it never fixes.
- `foam-scribe`: writes documentation, including the guide in section 8.

Step one is split into stages:

- **S1** Setup: fork detection, repository skeleton, `Allwmake`,
  `Alltest`, a (x, t) Gmsh case that meshes and imports.
- **S2** CCFV: `spaceTimeModel`, `cellCentred`, `analyticalSolution`,
  `spaceTimeAnalyticalFixedValue`, `spaceTimeErrors`, and the
  backward-Euler sanity check.
- **S3** VCFV geometry: `medianDualMesh`, the right-triangle unit test and
  the geometry identities, plus `perturbSpaceTimeMesh`.
- **S4** VCFV solver: `vertexCentred`, freestream and linear-exactness
  tests, and the Table 2 reproduction.
- **S5** Benchmark: the full refinement sweep, `.dat` tables and gnuplot
  figures.
- **S6** Guide: `docs/spaceTime4foamGuide.md`.

For every stage:

1. Present a short plan and wait for developer approval.
2. Delegate implementation to `foam-builder`.
3. Delegate review to `foam-auditor`.
4. If the verdict is FAIL, send the blocking issues back to
   `foam-builder`, then re-audit. Repeat until PASS.
5. Save the final audit report as `docs/audits/<stage>.md`.
6. Summarise the stage for the developer (what changed, test results,
   open questions) and wait before starting the next stage.

`foam-scribe` writes S6 once S5 data exist. Guide sections that do not
depend on results may be drafted earlier if the developer asks. The guide is
audited by `foam-auditor` like any other stage.

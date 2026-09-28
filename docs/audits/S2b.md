# S2b audit: spaceTimeLinearExtrapolation and S2 minors

Auditor: `foam-auditor` (read-only). Branch `step-one`, commits
`f0a0f2d..945c5a0` (8 commits, not pushed).

Environment: OpenFOAM v2412 (OpenFOAM.com), macOS, gmsh 4.15.2,
shellcheck 0.11.0, markdownlint-cli 0.49.1.

## Summary

**Verdict: PASS.**

- `spaceTimeLinearExtrapolation` implements the approved design and
  `CLAUDE.md` section 7: a `fixedGradient` condition with
  g = deltaCoeffs (grad(u)_P . (C_f - C_P)), so the face value uses the
  full offset vector. I checked this against the v2412 source myself.
- A clean rebuild has 0 compiler diagnostics. All 5 tests pass.
  shellcheck reports nothing on all 26 tracked scripts, and
  markdownlint reports nothing.
- The face-value test fails when I plant a normal-only fault or a
  sign fault.
- `ccfvOrder` reproduces the builder's tables exactly, wall time aside.
  CC-1, CC-2 and CC-2-LS are bitwise unchanged from S2.
- All four S2 minor items are fixed.
- I confirmed the builder's CC-2-LS-EX corner divergence, the gain of
  1.5, the Gauss gain of 1 and the N = 64 figures (2.6e33 with a
  residual of 8e-15). The top-right cell is the only cell with two EX
  faces.
- New finding: with Gauss linear the corner mode is exactly neutral.
  The discrete CC-2-EX system therefore has a one-parameter family of
  converged solutions.
  - With a = 1 this family changes only the two corner face values,
    never a cell value.
  - With a = 0.5 it changes the corner cell value and 109 cells near
    it.
  - From the standard `0/u` start the choice is deterministic (zero),
    which is why the benchmark numbers are reproducible.

  This is not blocking for the a = 1 benchmark, but it must be
  documented, and the developer must decide what to do before S5 (see
  Should fix 1 and 2).

## Blocking

None.

## Should fix

1. **The Gauss corner mode is neutral, so the solution is
   non-unique.** The header describes this incompletely.
   - Location:
     `spaceTimeLinearExtrapolationFvPatchScalarField.H`, lines 69-76.
   - Derivation (Gauss, any cell whose faces are all EX except one
     face D): the Gauss identity sum_f S_f d_f^T = V I gives the
     feedback F = I - S_D d_D^T / V, with d_D = C_D - C_P. Any
     gradient component v with d_D . v = 0 satisfies F v = v. In the
     top-right corner d_D is parallel to (1, 1), so the neutral
     direction is (1, -1). This holds for any a. It also holds on
     perturbed meshes, because the corner triangle has only boundary
     nodes, which `perturbSpaceTimeMesh` does not move.
   - Evidence, a = 1, N = 16. I restarted from the converged field
     with the corner face values changed by +0.01 (`tEnd`) and -0.01
     (`xRight`), which is pure (1, -1). The change is kept exactly for
     200 iterations. The residual stays at 5.2e-15 and no cell value
     changes (max |du| = 0).
   - Evidence, a = 0.5, N = 16. The same restart converges to a
     different solution with a residual of 4.6e-15. The corner cell
     value moves from 0.064579 to 0.062623, and 109 cells differ from
     the base run, by up to 1.2e-3 in the next cell (510).
   - The header says "with Gauss linear ... that component ... keeps
     its initial value (zero)". Its value is zero only when `0/u` has
     no `gradient` entry, because the zeroGradient start gives an
     exactly zero (1, -1) component. After a restart from a written
     state the component takes whatever value was written. "The cell
     balance does not see one direction" holds only for a = 1.
   - Consequence for the tables: the CC-2-EX corner face value is
     u_P + (the (1, 1) part only) . d. That gives the t = T error
     sin(pi h) at the corner (TExtrap Linf 0.02454 at N = 128, order
     1.000) and TExtrap L2 order 1.5.
   - Fix: document that the solution is non-unique for Gauss, the
     dependence on the start, and the a != 1 behaviour. Before S5 the
     developer should choose a corner treatment, for example
     zeroGradient on cells with two or more EX faces, or taking the
     corner-cell gradient from the upwind neighbour. Do not tune it;
     report the choice.
2. **CC-2-LS-EX is required by `CLAUDE.md` section 6 but is not
   delivered, and nothing stops a user from selecting it.**
   - The exclusion is justified (see the assessment below), but it
     changes the specification. It needs a developer decision and a
     `CLAUDE.md` update.
   - The tutorial `0/u` comment and the `fvSchemes` comment show how
     to switch EX on, and `fvSchemes` lists `leastSquares` for
     CC-2-LS. Neither says that EX with `leastSquares` is unsafe. A
     run with that combination reports "Converged" while the corner
     face values are wrong: at N = 16 check 2 differs by 0.100 at the
     1e-10 stop.
   - Fix: after the developer decides, either raise a
     `FatalError`/`Warning` in `updateCoeffs` when a cell has two or
     more faces with this condition and the gradient scheme is
     `leastSquares`, or at least add a warning comment in `0/u`.

## Minor

1. **Non-OpenFOAM.com branch.** The `#ifndef OPENFOAM_COM` branches use
   `getOrDefault`, `os.writeEntry` and `writeEntry("value", os)`. These
   are unguarded and untested on OpenFOAM.org and foam-extend. This is
   the same item as S2 Minor 3 ("open, accepted"). `IOobjectOption`
   and `writeValueEntry` also need a recent .com version.
2. **Check 1 uses the patch's own `gradSchemeName()`.** A wrong
   default name that happened to resolve to another existing scheme
   would not be caught. Low risk: I tested the option separately (see
   below).
3. **Check 2 (1e-8) is the only guard against a non-converged corner
   mode, and only `spaceTimeLinearExtrapolation` runs it (N = 16).**
   `ccfvOrder` has no face-value check for CC-2-EX. For Gauss this is
   acceptable, because I measured the distance from the 1e-10 stop to
   the fixed point (see the assessment below). Keep this in mind when
   perturbed meshes are added in S3 and S5.

## Scope and git

- The 8 commits are focused: the BC (`822222d`), the test (`16a2400`),
  `ccfvOrder` (`1aa0d19`), the tutorial (`141f383`) and one commit per
  minor (`da2983c`, `b8f5408`, `4b1de4c`, `945c5a0`). All carry the
  `Co-Authored-By` trailer.
- There is no VCFV, `medianDualMesh` or `perturbSpaceTimeMesh` code,
  and no change to `CLAUDE.md`, `docs/` or `.claude/`.
- `git status` was clean at the end of the audit.

## BC correctness

- **v2412 source** (`fixedGradientFvPatchField.C`):
  - `evaluate` (line 191) is
    `patchInternalField() + gradient_/patch().deltaCoeffs()`;
  - `valueInternalCoeffs` is 1 (line 205);
  - `valueBoundaryCoeffs` is `gradient()/deltaCoeffs()` (line 216).
- **deltaCoeffs** (`fvPatch.C:190`) is the boundary field of the
  mesh's `deltaCoeffs`, which is positive. With g = dc (grad . d), both
  the Gauss convection coefficients
  (`gaussConvectionScheme.C:107-108`) and `evaluate` give
  u_P + grad . d, up to the round-off of (dc x)/dc. The full vector d =
  `Cf() - Cn()` is used. `fvPatch::delta()` (line 160) would have kept
  only the normal part.
- **Update point.** The `fvMatrix` constructor calls
  `updateCoeffs()` (`fvMatrix.C:396`) before `fvmDiv` reads the
  coefficients. The gradient cache is keyed by name, and
  `gradScheme.C:122` checks it with `upToDate`.
- **Constructors.**
  - The dictionary constructor uses `LAZY_READ`. Without a `gradient`
    entry the face value is set to u_P and the gradient to zero; with
    one, the base class evaluates. Base-class construction dispatches
    to the base `updateCoeffs`, so this is safe.
  - The mapping, copy and iF constructors copy `gradSchemeName_`.
  - `autoMap` and `rmap` are inherited and map `gradient_`.
- **write.** It writes `gradient`, then `gradSchemeName` if it is not
  the default, then `value`. I checked this in a written file.
- **gradSchemeName option.** With `gradEX Gauss linear` and
  `gradSchemeName gradEX` on both patches, the run converges in 49
  iterations. The written field is bitwise identical to the test run,
  and the test passes.
- **Registration.** `TypeName("spaceTimeLinearExtrapolation")` and
  `makePatchTypeField` are present. `Make/files.openfoam` and
  `Make/files.foamextend` both list the source.
- **Cache claim.** With `cache` removed, CC-2-EX (N = 16) and CC-2
  (N = 64) are bitwise identical to the cached runs. Confirmed.
- **Lag.** Check 2 measures (grad u^{n+1} - grad u^n) . d. Its value
  of 6.45e-10 matches the per-iteration face change at the stop
  (1.1e-9 at iteration 49). The distance to the true fixed point is
  2.6e-9 on the faces at N = 16 and 8.3e-9 at N = 64. A tolerance of
  1e-8 is therefore justified at N = 16, but is tight at N = 64.
- **Planted faults** (scratch library builds, loaded first through
  `DYLD_LIBRARY_PATH`):

  | Fault | Check 1 | Check 2 | Exit |
  | --- | --- | --- | --- |
  | normal-only d (`patch().delta()`) | 5.36e-2 | 5.36e-2 | 1 |
  | wrong sign of g | 1.23e-1 | 1.23e-1 | 1 |

  A wrong C_P is unlikely to go unnoticed: the test uses
  `mesh.C()[faceCells]`, while the BC uses `patch().Cn()`.

## Assessment: is the CC-2-EX (Gauss) convergence claim sound?

It is sound for everything that iterates. The neutral component does
not iterate at all.

- I ran with tolerance 0 for 400 iterations, writing every iteration,
  and tracked u_P and both face values in the corner cell, and the
  maximum change of u:

  | N | stop | max du at stop | residual at 100 | max du at 100 |
  | --- | --- | --- | --- | --- |
  | 16 | 49 | 2.9e-10 (faces 1.1e-9) | 5.2e-15 | 0 (faces 1e-27) |
  | 64 | 46 | 1.1e-9 (faces 3.8e-9) | 6.3e-15 | 0 (faces 1e-29) |

- The iteration reaches an exact fixed point. The stop iterate is
  within 6.6e-10 (N = 16) and 2.3e-9 (N = 64) of it in the cells, and
  within 2.6e-9 and 8.3e-9 on the faces. That is negligible against
  errors of 3e-4 or more.
- The corner face values stay at u_P (about 1e-15) throughout. The
  exact values are -sin(pi h) on `tEnd` and +sin(pi h) on `xRight`.
- From a non-zero start the neutral mode does not drift: +/-0.01 is
  kept to all printed digits for 200 iterations (Should fix 1). It is
  neutral, not unstable.
- A perturbed internal field without a `gradient` entry cannot seed
  the mode. With u_f = u_P on both EX faces, the (1, -1) component of
  the Gauss gradient is exactly zero, and S_D . (1, -1) = 0.

## Assessment: CC-2-LS-EX corner mode

The builder's analysis is correct.

- **Geometry** (top-right cell on Left meshes): the vertices are
  (1-h, 1), (1, 1) and (1, 1-h), and C_P = (1-h/3, 1-h/3).
  - d_tEnd = (-h/6, h/3) and d_xRight = (h/3, -h/6).
  - The balance with a = 1 sees g . (d1 + d2) = g . (h/6, h/6) only.
- **Gain with leastSquares.** v2412 uses the normal-only offsets
  (0, h/3) and (h/3, 0), with weight 9/h^2, and the diagonal neighbour
  offset -(h/3)(1, 1), with weight 9/(2h^2).
  - The LS matrix is 0.5 [[1, 1], [1, 1]] + I, with eigenvalue 1 along
    (1, -1).
  - The EX face terms give (3/h)(g . d2, g . d1). Their (1, -1)
    component is (3/h)(h/2) g_w = 1.5 g_w.
  - The gain is therefore exactly 1.5, and it is (h/2)/(h/3) as the
    builder states.
- **Gain with Gauss:** exactly 1, from the identity in Should fix 1.
- **Reproduction** (tolerance 0, from the standard start):
  - N = 16: the mode is seeded by round-off and grows by 1.5 per
    iteration: 2.2e-14 at iteration 100, 1.8e-5 at 150, 1.1e4 at 200
    and 1.9e39 at 400. The residual falls below 1e-10 at iteration
    173 and returns to 0.199 by iteration 400.
  - N = 64: 2.5e-11 at 150, 1.6e-2 at 200 and 2.62e33 at 400, with a
    residual of 7.9e-15 at 400. This matches the builder.
  - A restart from the Gauss solution with +/-0.01 gives 0.0225,
    0.03375, 0.0506 and so on: a ratio of 1.5.
- **With the default 1e-10 stop:**
  - N = 16 converges in 173 iterations and check 2 fails with 0.100;
  - N = 64 converges in 152 iterations and all checks pass, because
    the stop comes before the round-off seed has grown. The stop hides
    the mode.
- **Other corners.** Parsing `owner` and `boundary` at N = 16 and 64,
  exactly two cells have two boundary faces:
  - cell 0: `tStart` and `xLeft`, both inflow;
  - the top-right cell: `tEnd` and `xRight`.

  The top-right cell is the only cell with two EX faces. On Right
  meshes no cell has two EX faces.
- **Open for S3/S5.** On perturbed meshes the corner triangle does not
  change, but the diagonal neighbour's centroid does. The LS gain then
  changes (not measured; `perturbSpaceTimeMesh` does not exist yet).
  The Gauss gain stays exactly 1.
- **Recommendation.** Excluding LS-EX is correct. It is not a tolerance
  problem: the fixed-point map is unstable in a direction that the
  residual cannot see. A fix needs a different corner treatment
  (Should fix 1), not relaxation.

## Reproducibility

- **Rerun against the builder.** After the clean rebuild, the
  `ccfvOrder` tables were compared with the builder's copies.
  `ccfvOrderOrders.dat` is identical. `ccfvOrder.dat` is identical
  apart from column 8 (wall time).
- **Against S2.** The CC-1, CC-2 and CC-2-LS rows and orders are
  identical to the S2 re-audit tables.
- **CC-2-EX claims** (all confirmed):
  - iterations 50, 49, 46, 46 and 61;
  - orders at 64 to 128: L1, L2 and Linf of 2.28, 2.25 and 2.21 (all
    cells), 2.28, 2.24 and 2.21 (interior) and 2.13, 2.16 and 2.24
    (boundary);
  - at N = 128, L2_all is 4.29e-4 against 1.71e-3 for CC-2 (4.0
    times lower) and Linf_all is 1.46e-3 against 2.94e-2 (20.1 times
    lower);
  - TExtrap Linf is 0.02454 = sin(pi/128), and TCell is order 1.02.

## S2 minor items

| Item | Evidence | Status |
| --- | --- | --- |
| `Allmesh 8 left extra` | usage message, exit 1 | Fixed |
| gmsh `diagonal foo` | "... not foo" error, exit 1 | Fixed |
| hasSource check | source 1e-9: new check FAILs | Fixed |
| blank lines | no run of 3 or more blank lines in any .C/.H | Fixed |

gmsh still writes `foo.msh` after the error, as the new `.geo`
comment states.

## Style

- All 19 tracked `.C` and `.H` files have the same 16-line licence
  header (one md5).
- The BC and the test follow OpenFOAM layout. The header documents the
  implementation, the lag, the cache and the corner limitation (but
  see Should fix 1).
- shellcheck (`.shellcheckrc` disables SC1091 only), on all 26 tracked
  scripts with a shebang (26 files have mode 100755): exit 0, no
  output.
- markdownlint on all tracked `.md` files: exit 0, no output.

## Commands run

```bash
git log/diff/show f0a0f2d..945c5a0; git status   # clean at end
./Allwclean; ./Allwmake     # exit 0; 0 lines with "warning:"/"error:"
./Alltest                   # 5 PASS, exit 0
diff ccfvOrder*.dat <builder copies>    # identical except wallTime
diff ccfvOrder*.dat <S2 copies>         # CC-1/2/2-LS identical
shellcheck <26 scripts>; markdownlint $(git ls-files '*.md')
# read in $FOAM_SRC: fixedGradientFvPatchField.C, fvPatch.C,
#   fvMatrix.C, gaussConvectionScheme.C, gradScheme.C
# scratch only (scratchpad/auditS2b):
#   CC-2-EX Gauss N = 16, 64: 400 iterations, tolerance 0, tracked
#   CC-2-LS-EX N = 16, 64: 400 iterations, tolerance 0, and at 1e-10
#   neutral-mode restarts (+/-0.01) for Gauss and LS, a = 1 and 0.5
#   planted faults: normal-only d, sign of g, source 1e-9
#   gradSchemeName option; cache removed; two-boundary-face cell count
#   Allmesh extra argument; gmsh -setstring diagonal foo
```

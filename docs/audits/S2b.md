# S2b audit: spaceTimeLinearExtrapolation, corner fallback, S2 minors

Auditor: `foam-auditor` (read-only). Branch `step-one`.

- **First audit:** `f0a0f2d..945c5a0`, which adds the BC, its test,
  CC-2-EX in `ccfvOrder` and the four S2 minor fixes.
- **Developer decision:** `da7606f` records the corner fallback in
  `CLAUDE.md` section 7.
- **Re-audit:** `da7606f..01fb434`, which adds the fallback, the
  guard, the extended test and CC-2-LS-EX in `ccfvOrder`.

Environment: OpenFOAM v2412 (OpenFOAM.com), macOS, gmsh 4.15.2,
shellcheck 0.11.0, markdownlint-cli 0.49.1. All runs were made in
`scratchpad/auditS2b/`. `git status` was clean at the end.

## Summary

**Verdict: PASS.**

- **BC.** `spaceTimeLinearExtrapolation` sets
  u_f = u_P + grad(u)_P . (C_f - C_P) with the full offset vector. It
  is a `fixedGradient` condition with g = deltaCoeffs (grad . d), and
  deltaCoeffs cancels exactly in both `evaluate` and the matrix
  coefficients (checked in the v2412 source).
- **Corner mode (first audit).** In the top-right cell, the only cell
  on `Left` meshes with two EX faces, the (1, -1) gradient component
  feeds back into itself:
  - with `leastSquares` the gain is exactly 1.5;
  - with Gauss the gain is exactly 1, a neutral mode, so the discrete
    solution is not unique.
- **Fallback (re-audit).** I show the mechanism directly in the tables
  below:
  - LS without the fallback: g_w grows by 1.500 per iteration while
    the residual stays near round-off;
  - LS with the fallback, from the same start: g_w stays at 1e-29;
  - a +/-0.01 corner perturbation is kept exactly without the fallback
    (Gauss) and removed in one update with it (Gauss and LS).
- **Rebuild and checks.** A clean rebuild has 0 compiler diagnostics
  and all 5 tests pass. shellcheck finds nothing in all 26 tracked
  scripts, and markdownlint finds nothing in all Markdown.
- **Reproducibility.**
  - CC-1, CC-2 and CC-2-LS are unchanged from S2.
  - CC-2-EX is unchanged from `945c5a0` in every error norm and in the
    iteration count. Only the final residual differs, in the 9th
    digit.
  - CC-2-LS-EX matches the builder's numbers exactly, wall time aside.
- **New finding (Should fix 1).** CC-2-LS-EX keeps a first-order
  outflow strip. Its Linf and L2 orders (1.00 and 1.56) come from the
  outflow boundary, not the inflow boundary. `leastSquares`
  extrapolation is not linearly exact there, so with LS the EX
  condition brings almost no gain over zeroGradient.
- **The 0.91 contraction of CC-2-LS-EX** is a local mode near the
  outflow corner. It is independent of N and of the fallback, and at
  N = 256 the run converges in 125 iterations. It is no risk for S5.

## Blocking

None.

## Should fix

1. **CC-2-LS-EX is first order in the outflow strip. The cause is
   `leastSquares` at the EX faces, not the inflow boundary.**
   - Evidence: max |e| by region at N = 128, from the `ccfvOrder`
     fields:

     | Scheme | tEnd strip | tStart strip | interior (> 3 h) |
     | --- | --- | --- | --- |
     | CC-2-LS-EX | 2.96e-2 | 2.56e-3 | 7.8e-4 |
     | CC-2-LS | 2.86e-2 | 2.56e-3 | 6.9e-4 |
     | CC-2-EX | 1.46e-3 | 4.7e-4 | 9.3e-4 |

     `xRight` gives the same values by the x <-> t symmetry. The
     CC-2-LS-EX tEnd-strip error is between 1.4e-2 and 3.0e-2 along
     the whole strip, not only at the corner.
   - Gradient error in the tEnd strip (N = 64, then N = 128):

     | Scheme | tEnd strip | interior |
     | --- | --- | --- |
     | CC-2-LS-EX | 13.2, then 13.4 (O(1)) | 0.26, then 0.22 |
     | CC-2-EX | 0.54, then 0.27 (O(h)) | 0.44, then 0.22 |

     For comparison, |grad u| = 8.9.
   - Mechanism: v2412 LS uses the normal-only boundary offset
     d_n = `fvPatch::delta()`, while EX sets u_f - u_P = g . d with the
     full d. For an exactly linear u, the LS fixed point is
     grad u + M^-1 w_f d_n (d - d_n) . grad u, so it is not exact. The
     tangential offset here is h/6.
   - The builder's attribution of the CC-2-LS-EX boundary order of
     about 1.0 to "the LS inflow inconsistency" is therefore incorrect
     for Linf and L2.
   - Fix: record this in `CLAUDE.md` section 7 and in the S6 guide.
     CC-2-LS-EX shows the S2 boundary-strip signature (L1 2.1, L2 1.5,
     Linf 1.0) because of the outflow faces. It is not a
     second-order outflow variant. Do not tune it.

## Minor

1. **The guard classifies schemes by name, not by type.**
   - Location: `leastSquaresGradient()` in the `.C` file. It treats a
     scheme as least-squares if any word contains "eastSquares" or
     equals "fourth".
   - I tested the stock v2412 scalar schemes with the fallback off and
     a +/-0.01 corner restart, using a no-guard build of `945c5a0`:

     | grad scheme | corner g_w over 30 iterations | guard |
     | --- | --- | --- |
     | `leastSquares` | x1.5 per iteration | Fatal |
     | `pointCellsLeastSquares` | converges to 6.70 | Fatal |
     | `edgeCellsLeastSquares` | converges to 6.70 | Fatal |
     | `iterativeGauss linear 5` | kept (-0.4526) | Warning |
     | `cellLimited Gauss linear 1` | kept (-0.4526) | Warning |

     The exact g_w at the corner is 2 pi sqrt(2) = 8.89.
   - For these schemes there is no false negative. The two stable LS
     variants are false positives: stopping them is conservative, but
     the message "gain above 1" is wrong for them.
   - A scheme from a user library whose name does not match gets only
     a Warning. `CLAUDE.md` says the guard stops "any configuration in
     which the unstable corner coupling could occur".
   - Recommendation: use an allow-list, so that only schemes whose
     first word is `Gauss` or `iterativeGauss`, or a limited Gauss
     scheme, give the Warning and all others are fatal. Alternatively,
     check the type of the constructed `fv::gradScheme`. The heuristic
     is acceptable in the meantime, because the fallback is on by
     default.
2. **Checks 5 and 6 do not test the BC's zeroing directly.**
   - Planted faults, run on the `spaceTimeLinearExtrapolation` test
     case (N = 16):

     | Fault | Gauss | leastSquares |
     | --- | --- | --- |
     | faces flagged, g not zeroed | all checks pass | checks 1, 2, 5 FAIL |
     | count per patch only | log 0; app passes | log 0; 1, 2, 5 FAIL |

   - The per-patch fault is caught by the solver-log count in
     `Allrun`, for both schemes.
   - The no-zeroing fault is caught only by LS. With Gauss and a = 1
     the corner gradient is zero by symmetry, so u_f = u_P anyway.
   - Check 6 counts the faces from the mesh, independently of the BC,
     so it tests the topology, not the BC. `-fallbackFaces 3` makes
     it FAIL (exit 1), so it can fail.
   - The coverage is adequate as long as the LS case stays in the test.
     Cheap hardening: also require the written `gradient` entry to be
     exactly 0 on the fallback faces. That discriminates for every
     scheme.
3. **Header wording.** The header says that without the fallback, for
   a != 1, "the corner cell value and cells downstream of it" depend on
   the start. The corner has no downstream cells. In my a = 0.5 test
   the 109 cells that changed lie upstream and near the corner (cells
   510, 509 and so on, decaying away from it), coupled through the
   gradients. Fix: say "cells near it".
4. **Non-OpenFOAM.com branch** (open from S2 Minor 3): `getOrDefault`,
   `os.writeEntry` and `writeEntry("value", os)` are unguarded and
   untested outside OpenFOAM.com.
5. **Check 1 uses the patch's own `gradSchemeName()`** (first audit).
   Low risk: the option was tested separately.

## First-audit items and their resolution

| First-audit item | Resolution |
| --- | --- |
| SF1 Gauss neutral mode, non-unique solution | Fallback; mechanism below |
| SF2 LS-EX required but missing, unguarded | Added with fallback; guard |
| M1 foam-extend/.org branch | Open (Minor 4) |
| M2 check 1 uses its own scheme name | Open (Minor 5) |
| M3 check 2 the only corner guard | Resolved by check 5 and LS case |
| S2 minors (Allmesh, .geo, hasSource, blanks) | Fixed (verified) |

Verification of the S2 minors:

- `Allmesh 8 left extra` prints a usage message and exits 1.
- gmsh with `-setstring diagonal foo` prints "diagonal must be left or
  right, not foo" and exits 1. It still writes the `.msh` file, as the
  `.geo` comment says.
- A planted `source()` of 1e-9 with `hasSource()` false makes the new
  check FAIL while "source against CD" still passes.
- No `.C` or `.H` file has a run of 3 or more blank lines.

## Mechanism: direct evidence

These runs are on the N = 16 `Left` mesh with a = 1, with the tolerance
set to 0 and every iteration written.

- The corner gradient is `postProcess -func "grad(u)"` on each written
  field. That is the gradient the BC uses in the next update.
- g_w = (g_x - g_y)/sqrt(2) is the (1, -1) projection.
- Runs without the fallback use a scratch library built from
  `945c5a0`, which has no guard. Runs with the fallback use the
  installed `01fb434` build.
- In every run below, the corner cell value u_P stays between 1e-17 and
  1e-14. The one exception is -2.9e-11 at iteration 200 without the
  fallback.

### leastSquares, standard start (`0/u`)

| k | g_w, no fallback | u_f tEnd, no fb | g_w, fallback | u_f tEnd, fb |
| --- | --- | --- | --- | --- |
| 1 | 0 | -4.9e-17 | 0 | -4.9e-17 |
| 50 | -3.0e-21 | -5.8e-15 | 8.1e-30 | -5.9e-15 |
| 100 | -1.90e-12 | 2.19e-14 | 7.3e-30 | -6.0e-15 |
| 120 | -6.315e-9 | 9.30e-11 | 7.3e-30 | -5.9e-15 |
| 121 | -9.473e-9 | 1.40e-10 | 7.5e-30 | -6.0e-15 |
| 150 | -1.2109e-3 | 1.78e-5 | 7.5e-30 | -5.9e-15 |
| 151 | -1.8164e-3 | 2.68e-5 | 7.5e-30 | -5.9e-15 |
| 175 | -30.58 | 0.450 | 7.3e-30 | -5.9e-15 |
| 200 | -7.72e5 | 1.14e4 | 7.3e-30 | -5.8e-15 |

- The ratios at 121/120 and 151/150 are 1.500. From 50 to 100 the
  average ratio is 1.50.
- Without the fallback, g_x = -g_y throughout, so the mode is pure
  (1, -1). u_f on `xRight` is -u_f on `tEnd`.
- The residual histories with and without the fallback are identical
  to the printed digits (4.715e-6 at k = 60 and 1.380e-8 at k = 120).
  The residual therefore does not see the mode at all.

### leastSquares, restart with a corner perturbation of +/-0.01

This starts from the converged state with the corner face values
changed by +0.01 on `tEnd` and -0.01 on `xRight`.

| k | g_x | g_y | g_w | u_f tEnd | u_f xRight |
| --- | --- | --- | --- | --- | --- |
| no fb 1 | -0.720 | 0.720 | -1.018 | 0.0150 | -0.0150 |
| no fb 2 | -1.080 | 1.080 | -1.527 | 0.0225 | -0.0225 |
| no fb 3 | -1.620 | 1.620 | -2.291 | 0.03375 | -0.03375 |
| no fb 5 | -3.645 | 3.645 | -5.155 | 0.0759 | -0.0759 |
| no fb 10 | -27.68 | 27.68 | -39.14 | 0.577 | -0.577 |
| no fb 40 | -5.31e6 | 5.31e6 | -7.51e6 | 1.11e5 | -1.11e5 |
| fb 1 | 4.2e-15 | 4.2e-15 | -1.6e-29 | -2.8e-15 | -2.8e-15 |
| fb 40 | -2.2e-15 | -2.2e-15 | 8.1e-30 | -5.9e-15 | -5.9e-15 |

Without the fallback the ratio is exactly 1.5 per iteration. With it,
the perturbation is removed in the first update.

### Gauss linear, restart with a corner perturbation of +/-0.01

| k | g_x | g_y | g_w | u_f tEnd | u_f xRight | residual |
| --- | --- | --- | --- | --- | --- | --- |
| no fb 1 | -0.320 | 0.320 | -0.4526 | 0.0100 | -0.0100 | 5.24e-15 |
| no fb 40 | -0.320 | 0.320 | -0.4526 | 0.0100 | -0.0100 | 5.24e-15 |
| fb 1 | 1.3e-15 | 1.3e-15 | 0 | -2.3e-15 | -2.3e-15 | 5.24e-15 |
| fb 40 | 1.3e-15 | 1.3e-15 | 0 | -2.3e-15 | -2.3e-15 | 5.24e-15 |

- Without the fallback the values are identical at every iteration
  from 1 to 40: the mode is neutral and invisible to the residual. The
  expected Gauss value is g_w = -4 (0.01)/(sqrt(2) h) = -0.4525.
- With the fallback the perturbation is removed in one update.
- At a = 0.5 without the fallback (first audit), the same restart
  converges to a different solution (residual 4.6e-15). The corner
  u_P moves from 0.064579 to 0.062623, and 109 cells change.

### Derivations (confirmed)

- **Gauss.** sum_f S_f d_f^T = V I, so the EX feedback is
  F = I - S_D d_D^T / V. Here d_D = -(h/6)(1, 1) is perpendicular to
  (1, -1), so F (1, -1) = (1, -1): the gain is 1 for any a and on
  perturbed meshes, because the corner triangle is fixed.
- **leastSquares.** The LS matrix is 0.5 [[1, 1], [1, 1]] + I, with
  eigenvalue 1 along (1, -1). The EX rows (9/h^2)(h/3) give 1.5 g_w,
  which is (h/2)/(h/3).
- **Other corners.** Parsing `owner` and `boundary` at N = 16 and 64,
  only two cells have two boundary faces: cell 0 (`tStart` and
  `xLeft`, both inflow) and the top-right cell. The fallback log and
  check 6 report 2 faces at every N, up to N = 256.

## Fallback implementation

- **Detection.** `multiFaces()` counts faces of this type per cell over
  all patches of the field, using `isA` and `faceCells`. It caches a
  `boolList` for the patch the first time it is called.
  - The cache is cleared in `autoMap` and `rmap`, and it is reset (not
    copied) in the copy and mapping constructors.
  - It is evaluated in `updateCoeffs` in both cases, fallback on and
    off, so the guard always runs.
- **Zeroing.** g is set to 0 on exactly the flagged faces, after the
  full-vector gradient is formed. Check 5 gives 0 difference from u_P
  on the fallback faces, and 15 of 16 faces per patch remain
  non-trivial (check 4).
- **Log.** There is one line per patch ("1 of N faces use the corner
  fallback") and exactly 2 lines per run, with no repeats from clones.
- **Switch.** `cornerFallback` is read with `getOrDefault<Switch>`
  (default true) and written only when it differs. In my runs:
  - `off` and `false` were read, written back as given, and gave the
    warning;
  - `yes` gave the fallback;
  - the default was not written.
- **Guard.** With the fallback off:
  - `leastSquares` gives a `FatalError` naming the patch and the
    scheme;
  - Gauss gives the neutral-mode Warning.

  The test really runs the solver for one iteration and requires a
  non-zero exit, "FOAM FATAL ERROR" and "is a least-squares scheme" in
  the log (confirmed in `run/CC-2-LS-EX-noFallback`).
- **Header and `CLAUDE.md`.** They agree on the rule, the reason
  (gain 1 and 1.5), the logging, the O(h) cost for a != 1 and the
  upgrade path (the upwind-neighbour gradient). The cost is visible:
  CC-2-EX at a = 0.5 has corner errors of 4.5e-3 (N = 64) and 2.7e-3
  (N = 128), and at N = 128 the corner cell sets Linf. See Minor 3 for
  the wording issue, and Minor 1 for the guard wording.

## BC correctness (first audit, still valid)

- **v2412 source** (`fixedGradientFvPatchField.C`):
  - `evaluate` (line 191) is `patchInternalField + gradient_/dc`;
  - `valueInternalCoeffs` is 1 (line 205);
  - `valueBoundaryCoeffs` is `gradient/dc` (line 216).
- **Offset.** The BC uses d = `Cf() - Cn()`. `fvPatch::delta()`
  (line 160) would have kept only the normal part.
- **Update point.** `fvMatrix.C:396` calls `updateCoeffs` before the
  convection coefficients are read (`gaussConvectionScheme.C:107`).
- **Read and write.** `LAZY_READ` means no `gradient` entry gives
  u_f = u_P. `gradSchemeName` is read and written correctly, and a
  named scheme gives a bitwise-identical field.
- **Cache claim.** With `cache` removed, CC-2-EX (N = 16) and CC-2
  (N = 64) are bitwise identical to the cached runs.
- **Lag and tolerance.** Check 2 gives 6.4e-10 for Gauss and 7.2e-10
  for LS. At the stop, the distance to the fixed point is 2.6e-9
  (N = 16) and 8.3e-9 (N = 64). The 1e-8 tolerance is justified at
  N = 16.
- **Planted faults in the BC** (first audit):
  - normal-only d gives check 1 = 5.36e-2, FAIL;
  - a sign flip gives check 1 = 1.23e-1, FAIL.

## Convergence claims

- **CC-2-EX (Gauss).** With the tolerance set to 0, the iteration
  reaches an exact fixed point: max |du| is 0 at iteration 100 for
  N = 16 and 64. At the 1e-10 stop the iterate is within 2.3e-9 of the
  fixed point in the cells. The residual criterion is therefore sound
  for everything that iterates, and the fallback removes the only
  component that did not iterate.
- **CC-2-LS-EX with the fallback.** The worst late-iteration residual
  ratio is 0.906, 0.910, 0.911, 0.911, 0.911 and 0.910 for
  N = 8, 16, 32, 64, 128 and 256.
  - The iterations fall with N: 174, 173, 164, 152, 139, and 125 at
    N = 256 (converged, loop wall time 0.84 s).
  - The slowest mode is local:
    - max |u^{k+1} - u^k| lies in the same cell relative to the corner
      at N = 16 and 64, at (0.67 h, 1.67 h), with a ratio of 0.90 in
      both;
    - its shape alternates in sign along the `tEnd` and `xRight`
      strips and decays within about 4 cells of the corner.
  - It is not caused by the fallback: the residual history without the
    fallback is identical.
  - It is benign. It is a mesh-scale mode with a fixed factor, so the
    work per decade of residual is constant in N. No risk at N = 256.
  - Compare CC-2-EX, whose worst ratio grows with N: 0.70 at N = 16
    and 0.82 at N = 128. This is the S2 deferred-correction behaviour.

## Reproducibility

- **Comparisons**, with the wall-time column removed:
  - The `ccfvOrder` rerun after the clean rebuild at `01fb434` matches
    my `945c5a0` rerun for CC-1, CC-2, CC-2-LS and CC-2-EX in every
    error norm and iteration count.
  - For CC-2-EX the final residual differs in the 9th digit
    (9.255515519e-11 against 9.255515582e-11 at N = 8), because the
    corner face values changed at round-off level.
  - `ccfvOrderOrders.dat` is identical for these schemes.
  - CC-1, CC-2 and CC-2-LS equal the S2 re-audit tables.
- **CC-2-LS-EX** equals the builder's log in every column except wall
  time.
  - Iterations: 174, 173, 164, 152 and 139.
  - Orders at 64 to 128: L1, L2 and Linf of 2.095, 1.559 and 1.005
    (all cells), and 1.009, 1.005 and 1.014 (boundary).

## Scope, git and style

- The 6 re-audit commits are focused and carry the trailer.
  `b046098` is the first-round report and `da7606f` is the developer's
  `CLAUDE.md` decision. There is no VCFV or S3 code.
- **Build.** `./Allwclean` and `./Allwmake` exit 0 with 0 lines
  matching `warning:` or `error:`.
- **Tests.** `./Alltest`: all 5 pass, exit 0. In
  `spaceTimeLinearExtrapolation`, checks 1-6 pass for CC-2-EX and
  CC-2-LS-EX, both guard checks pass, and the log counts are 2 and 2.
- **Licence header.** The same 16-line header in all 19 tracked `.C`
  and `.H` files (one md5).
- **Registration.** `TypeName` and `makePatchTypeField` are present,
  and both `Make/files*` lists include the BC.
- **shellcheck.** `.shellcheckrc` disables SC1091 only. On all 26
  tracked scripts with a shebang (all 26 have mode 100755): exit 0, no
  output.
- **markdownlint** on all tracked `.md` files: exit 0, no output.

## Recommendations to the developer

- Keep the fallback. It removes the only non-iterating or unstable
  component, and for a = 1 it changes no cell value.
- Report CC-2-LS-EX in S5 as an LS variant with a first-order outflow
  strip (Should fix 1), not as a fix for the outflow strip. CC-2-EX
  (Gauss) is the variant that removes the boundary-strip signature.
- Before any a != 1 case, apply the upgrade path: take the corner
  gradient from the upwind neighbour.

## Commands run

```bash
git log/diff/show f0a0f2d..945c5a0 and 945c5a0..01fb434; git status
./Allwclean; ./Allwmake            # twice (945c5a0, 01fb434): 0 diags
./Alltest                          # 5 PASS both times
diff ccfvOrder*.dat <945c5a0 rerun, S2 tables, builder log>
shellcheck <26 scripts>; markdownlint $(git ls-files '*.md')
# $FOAM_SRC read: fixedGradientFvPatchField.C, fvPatch.C, fvMatrix.C,
#   gaussConvectionScheme.C, gradScheme.C
# scratch only (scratchpad/auditS2b):
#   scratch lib from git archive 945c5a0 (no guard) for no-fallback runs
#   mech/: LS std start, LS and Gauss +/-0.01 restarts, fb and no fb,
#     grad(u) per iteration via postProcess; other grad schemes
#   fault libs: normal-only d, sign, source 1e-9 (first audit);
#     g not zeroed, per-patch count (re-audit); -fallbackFaces 3
#   error and gradient by region (N = 64, 128); a = 0.5 with fallback
#   LS-EX tolerance 0 at N = 16 and 64; LS-EX N = 256; switch values
```

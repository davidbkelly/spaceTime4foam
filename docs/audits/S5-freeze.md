# S5 regression freeze review

Reviewed on 6 October 2026 with OpenFOAM.com v2412 on macOS. Scope:
commits `4fbe3e8` through `66b3c77`, which implement S5 audit items
N1-N6 and the benchmark regression freeze. This review supplements
[S5.md](S5.md); it does not replace that numerical audit.

## Result

**PASS.** The committed reference was extracted by
`tests/benchmarkRegression/makeReference` from the accepted S5 results.
Default mode checks 214 N = 32 error values and two mesh checksums across
the two main families. Full mode adds 156 finest-pair orders and four
mesh checksums. CC-2-LS-EX on perturbed meshes remains a stability probe;
it has no accepted error or order value to freeze.

The reference uses relative tolerance `1e-4` for errors and absolute
tolerance `0.05` for orders. Values at the documented round-off level use
an absolute `1e-12` check. The default test is part of `./Alltest`; the
order check runs with `./Alltest -full benchmarkRegression`.

## Checks

- `./Alltest -full benchmarkRegression`: PASS, 216 default entries and
  160 full entries, with no failed comparisons.
- `./Alltest` after the environment correction: 20 passed, 0 failed.
- In a scratch copy of the regression runner, using the committed results
  and bypassing only the sweep, a 10% increase in the structured VC-2
  N = 32 `L2_all` value failed: `0.01030602135` expected versus
  `0.0113366` actual. A `+0.5` change in the perturbed VC-2
  128-to-256 `p_L2_all` value failed: `2.035` expected versus `2.535`
  actual. An unmodified control passed both modes. Scratch files and logs
  are outside the repository.
- `tests/gradientConsistency` checks the least-squares all-cell maxima at
  N = 16 and 64 on both families, covering N1. `tests/makeTables` plants
  `- - -`, `nan`, and swapped columns, and checks that a failed table
  generation leaves the old results intact, covering N2, N3 and N6.
- The figures use structured solid and perturbed dashed lines (N4).
  The probe watchdog computes its deadline with sub-second resolution
  (N5). `shellcheck` on the changed scripts and `markdownlint` on the
  existing documentation passed.
- A clean S5 re-sweep at commit `66b3c77` reproduced the accepted errors,
  observed orders, seed summaries, gradient diagnostic, probe outcomes,
  mesh checksums, and iteration counts exactly. Timing and provenance
  fields changed, as expected. The updated figures reflect the line-style
  change; the error and order data behind them did not change.

## Test environment correction

The first `./Alltest` run passed 19 tests and failed `makeTables`: macOS
created an `xcrun_db` cache under the test's `TMPDIR`, and the test counted
it as a leftover `makeTables` staging directory. The test now looks
specifically for `makeTables.*` directories. Its focused rerun passed.
The complete suite was then rerun with this correction and passed.

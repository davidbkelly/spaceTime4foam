# spaceTime4foam

`spaceTime4foam` is an OpenFOAM toolbox for finite volume methods on
space-time meshes. Its first benchmark solves one-dimensional travelling-wave
advection on a triangular `(x, t)` mesh and compares cell-centred (CCFV) and
vertex-centred (VCFV) discretisations. The same mesh is used for both methods.

The current implementation and numerical results have been tested with
**OpenFOAM.com v2412 in serial**. Time is a mesh coordinate; the solver's
`runTime` counts steady iterations rather than physical time.

## Documentation and results

- [Methods and results guide](docs/spaceTime4foamGuide.md)
  ([typeset PDF](docs/spaceTime4foamGuide.pdf))
- [Committed benchmark tables and figures](tutorials/advection1D/travellingSine/results/)
  and their [provenance manifest](tutorials/advection1D/travellingSine/results/manifest.dat)
- [Verification and audit reports](docs/audits/)
- [Reference citations](refs/README.md); the papers themselves are not
  included in the repository

## Build and verify

Install OpenFOAM.com v2412, Gmsh, and gnuplot. Source the OpenFOAM
environment, then run from the repository root:

```bash
source /path/to/OpenFOAM-v2412/etc/bashrc
./Allwmake
./Alltest
```

`Alltest` runs the verification suite and the frozen `N = 32` benchmark
checks. The optional full regression check also tests the finest-pair
observed orders:

```bash
./Alltest -full benchmarkRegression
```

## Run the benchmark

```bash
cd tutorials/advection1D/travellingSine
./Allrun
```

The full sweep covers structured and perturbed meshes and writes `.dat`
tables and PNG figures to `results/`. It can take substantial time. For a
single mesh, run `./Allmesh 32 left perturbed 12345` in the tutorial
directory. See the [guide](docs/spaceTime4foamGuide.md#7-code-map-and-how-to-run-it)
for scheme selection, a single-case workflow, and interpretation of the
results.

## Repository layout

| Path | Purpose |
| --- | --- |
| `applications/solvers/spaceTime4Foam/` | Steady space-time solver |
| `src/spaceTime4FoamModels/` | Models, geometry and boundary conditions |
| `applications/utilities/` | Mesh perturbation and error evaluation |
| `tutorials/advection1D/travellingSine/` | Benchmark and results |
| `tests/` | Verification and regression tests |
| `docs/` | Guide and independent audit reports |

The project is licensed under [GPL-3.0](LICENSE).

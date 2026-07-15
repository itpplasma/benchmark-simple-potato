# benchmark-simple-potato

This repository is the comparison and workflow layer for orbit benchmarks
between the sibling checkouts:

- `../SIMPLE`
- `../NEO-RT`

It is not the primary home for SIMPLE or NEO-RT implementation work. Use this
repo for:

- benchmark orchestration,
- shared Python environment setup,
- cross-code comparison scripts,
- benchmark inputs and run directories,
- plots, notes, and benchmark documentation.

## Repository purpose

The benchmark compares collisionless guiding-centre orbits in axisymmetric
tokamak fields:

- **SIMPLE**: symplectic, canonical-flux coordinates
- **POTATO / NEO-RT**: adaptive RK, cylindrical coordinates, finite-orbit-width
- **NEO-RT thin-orbit frequencies**: third reference where applicable

## Expected sibling layout

This repository assumes the following checkout structure:

```text
../benchmark-simple-potato
../SIMPLE
../NEO-RT
```

Many scripts and instructions assume those sibling paths exist.

## Shared Python environment

The recommended shared virtual environment for the whole benchmark workspace
lives in this repository:

```text
./.venv
```

Create it from the repository root with:

```bash
./setup-venv.sh
```

That script:

- creates or reuses `./.venv`,
- installs this repo's Python requirements,
- installs SIMPLE's Python requirements,
- installs `../SIMPLE` in editable mode,
- installs NEO-RT's Python requirements.

Reactivate later with:

```bash
source .venv/bin/activate
```

If the SIMPLE Python extension was previously built without Python support,
reactivate `.venv` first and reinstall SIMPLE from there.

### Portable builds on the shared filesystem

SIMPLE build directories under `/proj/plasma` are used from more than one
machine. Never build them with `CONFIG=Fast`: that profile enables
`-march=native`, so binaries or Python extensions can fail with an illegal
instruction when another machine accesses the same network directory. Use the
portable Release profile for shared builds:

```bash
source .venv/bin/activate
make -C ../SIMPLE clean
make -C ../SIMPLE CONFIG=Release
```

Use `Fast` only for a machine-local build directory that will never be executed
on another host.

## Where the important content lives

- `README.md` - entry point and current runnable workflow
- `WORKPLAN.md` - rung ladder, task split, status, and exit criteria
- `rung0/README.md` - equilibrium manifest, coordinate conventions, and exact
  Rung 0 commands
- `doc/benchmark.tex` - physics definitions, metrics, literature, and
  comparison specification
- `doc/refs.bib` - bibliography used by `doc/benchmark.tex`
- `rung0/` - public circular tokamak inputs and example run directories
- `doc/` - benchmark document source and generated LaTeX artifacts

The production Rung 0 equilibrium inputs are
`rung0/circ_chartmap_simple.nc` for SIMPLE and
`rung0/potato_run/circ.eqdsk` for POTATO. The chart-map radial coordinate is
`rho_tor = sqrt(s_tor)` and its geometry uses centimetres. SIMPLE commit
`08b85a1` or newer writes direct `R,Z` variables to `orbits.nc`; do not add a
canonical-to-VMEC conversion in benchmark tooling. Generated trajectories,
plots, and logs are run products, not equilibrium fixtures.

## Where to look next

### If the task is about the benchmark workflow

Start with:

1. `README.md`
2. `WORKPLAN.md`
3. `rung0/README.md`
4. `doc/benchmark.tex`
5. `doc/refs.bib` if citations or literature need updating

### If the task turns into a SIMPLE change

Switch to the sibling repo and read:

1. `../SIMPLE/CLAUDE.md`
2. `../SIMPLE/README.md`
3. `../SIMPLE/python/README.md`
4. `../SIMPLE/DOC/coordinates-and-fields.md`

### If the task turns into a NEO-RT change

Switch to the sibling repo and read:

1. `../NEO-RT/CLAUDE.md`
2. `../NEO-RT/README.md`
3. `../NEO-RT/doc/running.md`
4. `../NEO-RT/doc/file_formats.md`
5. `../NEO-RT/doc/library.md`

## Typical benchmark tasks

- set up or refresh the shared `.venv`,
- run or post-process SIMPLE and POTATO benchmark cases,
- compare trajectories, frequencies, or confinement outputs,
- update `WORKPLAN.md`,
- update `doc/benchmark.tex`,
- update `doc/refs.bib` when adding or changing citations,
- add benchmark-side tools or plotting scripts.

## Commands commonly relevant here

### Shared Python environment

```bash
./setup-venv.sh
source .venv/bin/activate
```

### Benchmark document

```bash
latexmk -pdf -cd doc/benchmark.tex
```

Run it from the repository root. Treat `benchmark.pdf`, `.aux`, `.log`, `.bbl`,
and `.fls` as build products.

## Working rules for this repo

- Treat this repo as the **integration / comparison** layer.
- Keep code changes for SIMPLE or NEO-RT in their own repositories unless the
  task is clearly benchmark-only.
- Keep notes, plots, and scripts here aligned with the current rung structure in
  `WORKPLAN.md`.
- Keep equilibrium provenance and run-directory contracts in
  `rung0/README.md`.
- When changing benchmark definitions, acceptance criteria, observables, or
  terminology, update both `WORKPLAN.md` and `doc/benchmark.tex` if they both
  describe the affected workflow.
- When documenting procedures, prefer sibling-relative paths (`../SIMPLE`,
  `../NEO-RT`) because that is the expected workspace layout.
- Keep the repository and `.venv` group read/writeable in the shared workspace.

## Maintenance

Update this file when any of the following change:

- sibling checkout assumptions,
- shared `.venv` behavior,
- main benchmark doc locations,
- benchmark workflow entrypoints,
- the relationship between this repo and the sibling SIMPLE / NEO-RT repos.
